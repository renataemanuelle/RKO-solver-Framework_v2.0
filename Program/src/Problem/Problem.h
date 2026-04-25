#include "readInstance.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <sstream>
#include <ctime>
#include <cmath>
#include <cctype>
#include <stdexcept>
#include <tuple>
#include <map>
#include <set>
#include <numeric>

//----------------- DEFINITION OF PROBLEM SPECIFIC TYPES -----------------------
struct TProblemData
{
    int n;

    int start_year;
    int start_month;
    int start_day;
    int start_hour;
    int start_minute;

    int horizon_hours;

    std::vector<int> norads;

    std::vector<Acquisition> acquisitions;

    // Pares stereo válidos (ângulo de visada entre 15° e 20°)
    std::vector<std::pair<int, int>> stereo_pairs;
    // Mapeamento: índice de aquisição → parceiros stereo válidos
    std::map<int, std::vector<int>> stereo_partners;

    // Métricas do cenário sobre dados brutos (pré-deduplicação),
    // para comparação direta com EOSPython (que não deduplica)
    int n_raw = 0;
    int requests_raw = 0;
    double sum_angle_raw = 0, sum_area_raw = 0, sum_price_raw = 0;
    double sum_sun_raw = 0, sum_cloud_raw = 0, sum_prio_raw = 0;
};


//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------

// Forward declarations
static void precompute_acquisition_fields(TProblemData& data, double satellite_height_km = 694.0);
static void compute_stereo_pairs(TProblemData& data);

/************************************************************************************
 Method: ReadData
 Description: read the input data
*************************************************************************************/
void ReadData(char name[], TProblemData &data)
{
    // Chama o leitor estruturado
    InstanceData instance = read_eos_instance(std::string(name));

    // Copiar dados básicos
    data.start_year = instance.start_year;
    data.start_month = instance.start_month;
    data.start_day = instance.start_day;
    data.start_hour = instance.start_hour;
    data.start_minute = instance.start_minute;

    data.horizon_hours = instance.horizon_hours;

    data.norads = instance.norads;

    data.acquisitions = instance.acquisitions;

    // Acumular métricas do cenário sobre dados brutos (pré-deduplicação)
    {
        std::set<std::string> raw_ids;
        for (const auto& a : data.acquisitions)
        {
            raw_ids.insert(a.ID);
            data.sum_angle_raw += a.angle;
            data.sum_area_raw  += a.area;
            data.sum_price_raw += a.price;
            data.sum_sun_raw   += a.sun_elevation;
            data.sum_cloud_raw += a.cloud_cover_real;
            data.sum_prio_raw  += a.priority;
        }
        data.n_raw = (int)data.acquisitions.size();
        data.requests_raw = (int)raw_ids.size();
    }

    // Deduplicação: remover linhas exatamente idênticas (artefato da expansão
    // stereo do LP no EOSPython — aquisições clonadas para simplificar S_constraint)
    {
        const int original_size = (int)data.acquisitions.size();
        std::set<std::string> seen;
        std::vector<Acquisition> unique;
        unique.reserve(original_size);

        for (auto& a : data.acquisitions)
        {
            std::string key = a.ID + "|" + std::to_string(a.satellite) + "|"
                            + a.time + "|" + a.satellite_location;
            if (seen.insert(key).second)
                unique.push_back(std::move(a));
        }

        data.acquisitions = std::move(unique);

        // Reatribuir index = posição no vetor (garante unicidade e rastreabilidade)
        for (int i = 0; i < (int)data.acquisitions.size(); i++)
            data.acquisitions[i].index = i;

        int removed = original_size - (int)data.acquisitions.size();
        if (removed > 0)
            std::cout << "Deduplication: removed " << removed
                      << " duplicate rows (" << original_size << " -> "
                      << data.acquisitions.size() << ")\n";
    }

    // Tamanho do vetor random-key
    data.n = data.acquisitions.size();

    // Impressão para confirmação
    std::cout << "\n=== EOS DATA LOADED (via readInstance) ===\n";
    std::cout << "Start: " << data.start_year << "-" << data.start_month << "-" << data.start_day
          << " " << data.start_hour << ":" << data.start_minute << "\n";
    std::cout << "Horizon (hours): " << data.horizon_hours << "\n";
    std::cout << "NORAD count: " << data.norads.size() << "\n\n";
    std::cout << "Total acquisitions: " << data.n << "\n";

    // Pré-calcular campos derivados (epoch_seconds, sat_xyz, req_xyz)
    precompute_acquisition_fields(data);

    // Pré-computar pares stereo válidos
    compute_stereo_pairs(data);
    std::cout << "Stereo pairs found: " << data.stereo_pairs.size() << "\n";
}


//-------------------------- AUXILIARY TYPES FOR PROBLEM SOLUTION --------------------------
struct DecodedSolution
{
    std::vector<int> x;                 // 0/1 por aquisição
    std::vector<int> selected_idxs;     // índices selecionados
    double objective_value;             // por enquanto, soma dos score_scenario
};

struct Interval
{
    long long start; // epoch seconds
    long long end;   // epoch seconds
};

static constexpr double PI = 3.14159265358979323846;

//-------------------------- AUXILIARY FUNCTIONS FOR PROBLEM SOLUTION --------------------------

// Converte graus para radianos
static inline double deg2rad(double deg)
{
    return deg * PI / 180.0;
}

// Remove "np.float64(" e ")" da string, se existirem
static std::string clean_location_string(std::string s)
{
    const std::string token = "np.float64(";

    while (true)
    {
        size_t pos = s.find(token);
        if (pos == std::string::npos) break;
        s.erase(pos, token.size());

        size_t close = s.find(')', pos);
        if (close != std::string::npos) s.erase(close, 1);
    }

    return s;
}

// Extrai latitude e longitude de strings como:
// "[57.50037782967006, 9.785145863466415]"
// "[np.float64(55.64321798038925), np.float64(12.16354832122793)]"
static std::pair<double, double> parse_lat_lon(const std::string& loc)
{
    std::string s = clean_location_string(loc);

    // manter apenas números, sinais, ponto, vírgula e expoente
    for (char& c : s)
    {
        if (!(std::isdigit(static_cast<unsigned char>(c)) ||
              c == '-' || c == '+' || c == '.' || c == ',' ||
              c == 'e' || c == 'E'))
        {
            c = ' ';
        }
    }

    std::replace(s.begin(), s.end(), ',', ' ');

    std::stringstream ss(s);
    double lat, lon;
    ss >> lat >> lon;

    if (ss.fail())
        throw std::runtime_error("Erro ao converter latitude/longitude: " + loc);

    return {lat, lon};
}

// Converte latitude/longitude para cartesiano
// elevation_km = 0 para request, = altura do satélite para satélite
static Vec3 cart_system(double lat_deg, double lon_deg, double elevation_km)
{
    const double R = 6371.0 + elevation_km;
    const double lat = deg2rad(lat_deg);
    const double lon = deg2rad(lon_deg);

    return {
        R * std::cos(lat) * std::cos(lon),
        R * std::cos(lat) * std::sin(lon),
        R * std::sin(lat)
    };
}

// Subtrai dois vetores
static inline Vec3 subtract(const Vec3& a, const Vec3& b)
{
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

// Calcula o produto escalar de dois vetores
static inline double dot(const Vec3& a, const Vec3& b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Calcula a norma de um vetor
static inline double norm(const Vec3& v)
{
    return std::sqrt(dot(v, v));
}

// Converte uma string de tempo para segundos desde a época (epoch)
static long long parse_time_to_epoch_seconds(const std::string& time_str)
{
    std::tm tm = {};
    std::istringstream ss(time_str);

    // exemplo esperado: "2024-01-01 12:34:56"
    ss >> std::get_time(&tm, "%Y-%m-%d %H:%M:%S");
    if (ss.fail())
        throw std::runtime_error("Erro ao converter tempo: " + time_str);

    // mktime usa horário local; se você já tiver outra função no projeto, use a sua.
    return static_cast<long long>(std::mktime(&tm));
}

// Calcula o ângulo entre as duas linhas de visada, em graus.
// Usa coordenadas cartesianas pré-calculadas (sat_xyz, req_xyz).
static double maneuver_angle_deg_fast(const Acquisition& a, const Acquisition& b)
{
    const Vec3 vec1 = subtract(a.sat_xyz, a.req_xyz);
    const Vec3 vec2 = subtract(b.sat_xyz, b.req_xyz);

    const double n1 = norm(vec1);
    const double n2 = norm(vec2);

    if (n1 == 0.0 || n2 == 0.0)
        throw std::runtime_error("Norma zero ao calcular vetor de observação.");

    double cos_theta = dot(vec1, vec2) / (n1 * n2);
    cos_theta = std::clamp(cos_theta, -1.0, 1.0);

    return std::acos(cos_theta) * 180.0 / PI;
}

// Verifica se a sequência a -> b é viável por manobrabilidade.
// Usa epoch_seconds e coordenadas pré-calculadas.
static bool maneuver_feasible(const Acquisition& a,
                              const Acquisition& b,
                              double rotation_speed_deg_per_sec = 30.0 / 12.0)
{
    if (a.satellite != b.satellite)
        return false;

    const double delta_t = static_cast<double>(b.epoch_seconds - a.epoch_seconds);
    if (delta_t < 0.0)
        return false;

    const double angle_deg = maneuver_angle_deg_fast(a, b);
    const double t_man = angle_deg / rotation_speed_deg_per_sec;

    return (a.duration + t_man) <= delta_t;
}

// Pré-calcula epoch_seconds e coordenadas cartesianas para todas as aquisições.
// Deve ser chamada uma vez após a deduplicação, antes de compute_stereo_pairs.
static void precompute_acquisition_fields(TProblemData& data,
                                          double satellite_height_km)
{
    for (auto& a : data.acquisitions)
    {
        std::tm tm = {};
        std::istringstream ss(a.time);
        ss >> std::get_time(&tm, "%Y-%m-%d %H:%M:%S");
        a.epoch_seconds = static_cast<long long>(std::mktime(&tm));

        auto [slat, slon] = parse_lat_lon(a.satellite_location);
        a.sat_xyz = cart_system(slat, slon, satellite_height_km);

        auto [rlat, rlon] = parse_lat_lon(a.request_location);
        a.req_xyz = cart_system(rlat, rlon, 0.0);
    }
}

// Constantes para validação de pares stereo (conforme EOSPython)
static constexpr double STEREO_ANGLE_CENTER = 17.5;
static constexpr double STEREO_ANGLE_ERROR  = 2.5;

// Pré-computa pares stereo válidos: para cada ID com stereo > 0,
// verifica se o ângulo entre as linhas de visada está em [15°, 20°]
static void compute_stereo_pairs(TProblemData& data)
{
    data.stereo_pairs.clear();
    data.stereo_partners.clear();

    // Agrupar aquisições stereo por ID
    std::map<std::string, std::vector<int>> stereo_groups;
    for (int i = 0; i < data.n; i++)
    {
        if (data.acquisitions[i].stereo > 0)
            stereo_groups[data.acquisitions[i].ID].push_back(i);
    }

    // Para cada grupo, verificar todos os pares
    for (const auto& [id, indices] : stereo_groups)
    {
        if (indices.size() <= 1) continue;

        for (size_t a = 0; a < indices.size() - 1; a++)
        {
            for (size_t b = a + 1; b < indices.size(); b++)
            {
                double angle = maneuver_angle_deg_fast(
                    data.acquisitions[indices[a]],
                    data.acquisitions[indices[b]]);

                if (angle >= STEREO_ANGLE_CENTER - STEREO_ANGLE_ERROR &&
                    angle <= STEREO_ANGLE_CENTER + STEREO_ANGLE_ERROR)
                {
                    data.stereo_pairs.push_back({(int)indices[a], (int)indices[b]});
                    data.stereo_partners[(int)indices[a]].push_back((int)indices[b]);
                    data.stereo_partners[(int)indices[b]].push_back((int)indices[a]);
                }
            }
        }
    }
}

static bool can_insert_by_maneuver(
    const std::vector<int>& selected_idxs,
    int cand_idx,
    const TProblemData& data)
{
    const Acquisition& cand = data.acquisitions[cand_idx];
    const long long cand_time = cand.epoch_seconds;

    // Busca binária pela posição temporal de inserção — O(log n)
    auto it = std::lower_bound(
        selected_idxs.begin(), selected_idxs.end(), cand_time,
        [&](int idx, long long t) {
            return data.acquisitions[idx].epoch_seconds < t;
        });
    int pos = (int)(it - selected_idxs.begin());

    if (pos > 0)
    {
        if (!maneuver_feasible(data.acquisitions[selected_idxs[pos - 1]], cand))
            return false;
    }

    if (pos < (int)selected_idxs.size())
    {
        if (!maneuver_feasible(cand, data.acquisitions[selected_idxs[pos]]))
            return false;
    }

    return true;
}

// Insere uma aquisição na lista ordenada por tempo — O(log n) para busca
static void insert_sorted_by_time(
    std::vector<int>& selected_idxs,
    int cand_idx,
    const TProblemData& data)
{
    const long long cand_time = data.acquisitions[cand_idx].epoch_seconds;

    auto it = std::lower_bound(
        selected_idxs.begin(), selected_idxs.end(), cand_time,
        [&](int idx, long long t) {
            return data.acquisitions[idx].epoch_seconds < t;
        });

    selected_idxs.insert(it, cand_idx);
}

// Verifica se um intervalo pode ser inserido sem sobreposição temporal (não modifica o vetor)
static bool can_insert_interval_no_overlap(const std::vector<Interval>& used, const Interval& cand)
{
    auto it = std::lower_bound(
        used.begin(), used.end(), cand.start,
        [](const Interval& a, long long s){ return a.start < s; }
    );

    if (it != used.begin())
    {
        const Interval& prev = *(it - 1);
        if (cand.start < prev.end) return false;
    }

    if (it != used.end())
    {
        const Interval& next = *it;
        if (cand.end > next.start) return false;
    }

    return true;
}

// Insere intervalo mantendo ordenação por start (chamar após verificar com can_insert_interval_no_overlap)
static void insert_interval_sorted(std::vector<Interval>& used, const Interval& cand)
{
    auto it = std::lower_bound(
        used.begin(), used.end(), cand.start,
        [](const Interval& a, long long s){ return a.start < s; }
    );
    used.insert(it, cand);
}

// Decodifica uma solução random-key para uma solução do problema
static DecodedSolution decode_solution(
    const TSol& s,
    const TProblemData& data)
{
    DecodedSolution out;
    out.x.assign(data.n, 0);
    out.objective_value = 0.0;

    std::map<int, std::vector<Interval>> sat_intervals;
    std::map<std::string, int> id_selection_count;

    // Tenta inserir uma aquisição verificando todas as restrições.
    // Retorna true se já estava selecionada ou se foi inserida com sucesso.
    auto try_insert = [&](int idx) -> bool {
        if (out.x[idx] == 1) return true;

        const Acquisition& acq = data.acquisitions[idx];

        int max_per_id = std::max(acq.stereo + 1, acq.strips);
        if (id_selection_count[acq.ID] >= max_per_id)
            return false;

        Interval interval = {acq.epoch_seconds,
                              acq.epoch_seconds + static_cast<long long>(std::ceil(acq.duration))};

        if (!can_insert_interval_no_overlap(sat_intervals[acq.satellite], interval))
            return false;

        if (!can_insert_by_maneuver(out.selected_idxs, idx, data))
            return false;

        id_selection_count[acq.ID]++;
        insert_interval_sorted(sat_intervals[acq.satellite], interval);
        insert_sorted_by_time(out.selected_idxs, idx, data);
        out.x[idx] = 1;
        out.objective_value += acq.score_scenario;
        return true;
    };

    // Desfaz inserção de uma aquisição (inversa de try_insert)
    auto undo_insert = [&](int idx) {
        const Acquisition& acq = data.acquisitions[idx];
        out.x[idx] = 0;
        out.objective_value -= acq.score_scenario;
        id_selection_count[acq.ID]--;

        auto sit = std::find(out.selected_idxs.begin(), out.selected_idxs.end(), idx);
        if (sit != out.selected_idxs.end()) out.selected_idxs.erase(sit);

        auto& intervals = sat_intervals[acq.satellite];
        auto iit = std::find_if(intervals.begin(), intervals.end(),
            [&acq](const Interval& iv) { return iv.start == acq.epoch_seconds; });
        if (iit != intervals.end()) intervals.erase(iit);
    };

    std::vector<int> idx(data.n);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(),
              [&](int a, int b)
              {
                  return s.rk[a] > s.rk[b];
              });

    // Inserção gulosa com enforcement eager de pares stereo (conforme EOSPython)
    for (int k = 0; k < data.n; k++)
    {
        int cand_idx = idx[k];

        if (out.x[cand_idx] == 1) continue;

        const Acquisition& cand = data.acquisitions[cand_idx];

        auto partners_it = data.stereo_partners.find(cand_idx);
        bool has_valid_partners = (cand.stereo > 0) &&
            (partners_it != data.stereo_partners.end()) &&
            !partners_it->second.empty();

        if (has_valid_partners)
        {
            bool partner_already_selected = false;
            for (int p : partners_it->second)
            {
                if (out.x[p] == 1) { partner_already_selected = true; break; }
            }

            if (partner_already_selected)
            {
                try_insert(cand_idx);
            }
            else
            {
                if (!try_insert(cand_idx))
                    continue;

                std::vector<int> partners = partners_it->second;
                std::sort(partners.begin(), partners.end(), [&](int a, int b) {
                    return data.acquisitions[a].score_scenario
                         > data.acquisitions[b].score_scenario;
                });

                bool found_partner = false;
                for (int p : partners)
                {
                    if (out.x[p] == 1) { found_partner = true; break; }
                    if (try_insert(p))  { found_partner = true; break; }
                }

                if (!found_partner)
                    undo_insert(cand_idx);
            }
        }
        else
        {
            try_insert(cand_idx);
        }
    }

    return out;
}

/************************************************************************************
 Method: Decoder 
 Description: mapping the random-key solution into a problem solution
*************************************************************************************/
double Decoder(TSol &s, const TProblemData &data)
{
    DecodedSolution decoded = decode_solution(s, data);

    s.selected_idxs = decoded.selected_idxs;
    s.x = decoded.x;

    // Negamos o score porque o framework RKO MINIMIZA, mas o EOS quer MAXIMIZAR.
    // Assim, score maior → ofv menor → melhor no ranking do pool.
    return -decoded.objective_value;
}


/************************************************************************************
 Method: FreeMemoryProblem
 Description: Free local memory allocate by Problem
*************************************************************************************/
void FreeMemoryProblem(TProblemData &data){
    data.norads.clear();
    data.acquisitions.clear();
}
