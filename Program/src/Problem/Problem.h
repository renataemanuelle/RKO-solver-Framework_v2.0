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
};


//-------------------------- FUNCTIONS OF SPECIFIC PROBLEM --------------------------

// Forward declaration (definida após maneuver_angle_deg)
static void compute_stereo_pairs(TProblemData& data, double satellite_height_km = 694.0);

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

    // Tamanho do vetor random-key
    data.n = data.acquisitions.size();

    // Impressão para confirmação
    std::cout << "\n=== EOS DATA LOADED (via readInstance) ===\n";
    std::cout << "Start: " << data.start_year << "-" << data.start_month << "-" << data.start_day
          << " " << data.start_hour << ":" << data.start_minute << "\n";
    std::cout << "Horizon (hours): " << data.horizon_hours << "\n";
    std::cout << "NORAD count: " << data.norads.size() << "\n\n";
    std::cout << "Total acquisitions: " << data.n << "\n";

    int to_print = std::min(3, data.n);

    for (int i = 0; i < to_print; i++)
    {
        const Acquisition& a = data.acquisitions[i];

        std::cout << "\n--- Acquisition " << i << " ---\n";
        std::cout << "index: " << a.index << "\n";
        std::cout << "ID: " << a.ID << "\n";
        std::cout << "stereo: " << a.stereo << "\n";
        std::cout << "satellite: " << a.satellite << "\n";
        std::cout << "satellite_location: " << a.satellite_location << "\n";
        std::cout << "request_location: " << a.request_location << "\n";
        std::cout << "time: " << a.time << "\n";
        std::cout << "area: " << a.area << "\n";
        std::cout << "strips: " << a.strips << "\n";
        std::cout << "duration: " << a.duration << "\n";
        std::cout << "distance: " << a.distance << "\n";
        std::cout << "angle: " << a.angle << "\n";
        std::cout << "sun_elevation: " << a.sun_elevation << "\n";
        std::cout << "cloud_cover_estimate: " << a.cloud_cover_estimate << "\n";
        std::cout << "priority: " << a.priority << "\n";
        std::cout << "priority_mod: " << a.priority_mod << "\n";
        std::cout << "customer_type_mod: " << a.customer_type_mod << "\n";
        std::cout << "price: " << a.price << "\n";
        std::cout << "waiting_time: " << a.waiting_time << "\n";
        std::cout << "uncertainty: " << a.uncertainty << "\n";
        std::cout << "cloud_cover_real: " << a.cloud_cover_real << "\n";
        std::cout << std::fixed << std::setprecision(16);
        std::cout << "score_scenario: " << a.score_scenario << "\n";
        std::cout << "score_method: " << a.score_method << "\n";
        std::cout << "score_alpha: " << a.score_alpha << "\n";
    }

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

struct Vec3
{
    double x, y, z;
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

// Calcula o ângulo entre as duas linhas de visada, em graus
static double maneuver_angle_deg(const Acquisition& a,
                                 const Acquisition& b,
                                 double satellite_height_km)
{
    auto [sat_lat1, sat_lon1] = parse_lat_lon(a.satellite_location);
    auto [req_lat1, req_lon1] = parse_lat_lon(a.request_location);

    auto [sat_lat2, sat_lon2] = parse_lat_lon(b.satellite_location);
    auto [req_lat2, req_lon2] = parse_lat_lon(b.request_location);

    const Vec3 sat_xyz1 = cart_system(sat_lat1, sat_lon1, satellite_height_km);
    const Vec3 req_xyz1 = cart_system(req_lat1, req_lon1, 0.0);
    const Vec3 vec1 = subtract(sat_xyz1, req_xyz1);

    const Vec3 sat_xyz2 = cart_system(sat_lat2, sat_lon2, satellite_height_km);
    const Vec3 req_xyz2 = cart_system(req_lat2, req_lon2, 0.0);
    const Vec3 vec2 = subtract(sat_xyz2, req_xyz2);

    const double n1 = norm(vec1);
    const double n2 = norm(vec2);

    if (n1 == 0.0 || n2 == 0.0)
        throw std::runtime_error("Norma zero ao calcular vetor de observação.");

    double cos_theta = dot(vec1, vec2) / (n1 * n2);
    cos_theta = std::clamp(cos_theta, -1.0, 1.0);

    return std::acos(cos_theta) * 180.0 / PI;
}

// MÉTODO PRINCIPAL:
// verifica se a sequência a -> b é viável por manobrabilidade
static bool maneuver_feasible(const Acquisition& a,
                              const Acquisition& b,
                              double satellite_height_km = 694.0,
                              double rotation_speed_deg_per_sec = 30.0 / 12.0)
{
    // no EOSPython a checagem é feita entre tentativas do mesmo satélite
    if (a.satellite != b.satellite)
        return false;

    const long long ta = parse_time_to_epoch_seconds(a.time);
    const long long tb = parse_time_to_epoch_seconds(b.time);

    const double delta_t = static_cast<double>(tb - ta);
    if (delta_t < 0.0)
        return false;

    const double angle_deg = maneuver_angle_deg(a, b, satellite_height_km);
    const double t_man = angle_deg / rotation_speed_deg_per_sec;

    return (a.duration + t_man) <= delta_t;
}

// Constantes para validação de pares stereo (conforme EOSPython)
static constexpr double STEREO_ANGLE_CENTER = 17.5;
static constexpr double STEREO_ANGLE_ERROR  = 2.5;

// Pré-computa pares stereo válidos: para cada ID com stereo > 0,
// verifica se o ângulo entre as linhas de visada está em [15°, 20°]
static void compute_stereo_pairs(TProblemData& data, double satellite_height_km)
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
                double angle = maneuver_angle_deg(
                    data.acquisitions[indices[a]],
                    data.acquisitions[indices[b]],
                    satellite_height_km);

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
    const long long cand_time = parse_time_to_epoch_seconds(cand.time);

    // encontrar posição temporal de inserção
    int pos = 0;
    while (pos < (int)selected_idxs.size())
    {
        const Acquisition& cur = data.acquisitions[selected_idxs[pos]];
        long long cur_time = parse_time_to_epoch_seconds(cur.time);

        if (cand_time < cur_time)
            break;

        pos++;
    }

    // checa com predecessor
    if (pos > 0)
    {
        const Acquisition& prev = data.acquisitions[selected_idxs[pos - 1]];
        if (!maneuver_feasible(prev, cand))
            return false;
    }

    // checa com sucessor
    if (pos < (int)selected_idxs.size())
    {
        const Acquisition& next = data.acquisitions[selected_idxs[pos]];
        if (!maneuver_feasible(cand, next))
            return false;
    }

    return true;
}

// Insere uma aquisição na lista ordenada por tempo
static void insert_sorted_by_time(
    std::vector<int>& selected_idxs,
    int cand_idx,
    const TProblemData& data)
{
    const long long cand_time =
        parse_time_to_epoch_seconds(data.acquisitions[cand_idx].time);

    int pos = 0;
    while (pos < (int)selected_idxs.size())
    {
        const long long cur_time =
            parse_time_to_epoch_seconds(data.acquisitions[selected_idxs[pos]].time);

        if (cand_time < cur_time)
            break;

        pos++;
    }

    selected_idxs.insert(selected_idxs.begin() + pos, cand_idx);
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

    // Intervalos ocupados por satélite (para detectar sobreposição temporal)
    std::map<int, std::vector<Interval>> sat_intervals;

    // Contador de seleções por ID: limite = max(stereo+1, strips) conforme EOSPython
    std::map<std::string, int> id_selection_count;

    // 1) criar lista de índices
    std::vector<int> idx(data.n);
    std::iota(idx.begin(), idx.end(), 0);

    // 2) ordenar por random-key decrescente (maior chave = maior prioridade)
    std::sort(idx.begin(), idx.end(),
              [&](int a, int b)
              {
                  return s.rk[a] > s.rk[b];
              });

    // 3) inserção gulosa com verificação de restrições
    for (int k = 0; k < data.n; k++)
    {
        int cand_idx = idx[k];
        const Acquisition& cand = data.acquisitions[cand_idx];

        // Restrição 1: limite de seleções por ID
        // max(stereo+1, strips): normal→1, stereo=1→2, strips=2→2, strips=3→3
        int max_per_id = std::max(cand.stereo + 1, cand.strips);
        if (id_selection_count[cand.ID] >= max_per_id)
            continue;

        // Intervalo temporal da aquisição candidata
        long long cand_start = parse_time_to_epoch_seconds(cand.time);
        long long cand_end = cand_start + static_cast<long long>(std::ceil(cand.duration));
        Interval cand_interval = {cand_start, cand_end};

        // Restrição 2: sobreposição temporal no mesmo satélite
        if (!can_insert_interval_no_overlap(sat_intervals[cand.satellite], cand_interval))
            continue;

        // Restrição 3: viabilidade de manobra com vizinhos temporais
        if (!can_insert_by_maneuver(out.selected_idxs, cand_idx, data))
            continue;

        // Todas as restrições satisfeitas: inserir a aquisição
        id_selection_count[cand.ID]++;
        insert_interval_sorted(sat_intervals[cand.satellite], cand_interval);
        insert_sorted_by_time(out.selected_idxs, cand_idx, data);
        out.x[cand_idx] = 1;
        out.objective_value += cand.score_scenario;
    }

    // Pós-processamento: restrição de pares stereo (ambos selecionados ou nenhum)
    // Para cada par válido, se apenas um foi selecionado, remover o selecionado.
    // Iterar até estabilizar (remoção de um pode afetar pares encadeados).
    bool any_stereo_removal = false;
    {
        bool changed = true;
        while (changed)
        {
            changed = false;
            for (const auto& [idx_a, idx_b] : data.stereo_pairs)
            {
                if (out.x[idx_a] != out.x[idx_b])
                {
                    int to_remove = out.x[idx_a] ? idx_a : idx_b;
                    out.x[to_remove] = 0;
                    out.objective_value -= data.acquisitions[to_remove].score_scenario;
                    changed = true;
                    any_stereo_removal = true;
                }
            }
        }
    }

    // Se houve remoções no pós-processamento, reconstruir selected_idxs
    if (any_stereo_removal)
    {
        out.selected_idxs.clear();
        for (int i = 0; i < data.n; i++)
        {
            if (out.x[i]) out.selected_idxs.push_back(i);
        }
        std::sort(out.selected_idxs.begin(), out.selected_idxs.end(),
            [&](int a, int b) {
                return parse_time_to_epoch_seconds(data.acquisitions[a].time)
                     < parse_time_to_epoch_seconds(data.acquisitions[b].time);
            });
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
