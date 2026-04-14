#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <iostream>

// -----------------------------
// Tipos básicos do leitor
// -----------------------------
struct DateTime {
    int year = 0, month = 0, day = 0, hour = 0, minute = 0;

    bool valid() const {
        return year > 0 && month >= 1 && month <= 12 && day >= 1 && day <= 31 &&
               hour >= 0 && hour <= 23 && minute >= 0 && minute <= 59;
    }
};

struct Vec3 { double x, y, z; };

struct Acquisition
{
    int index;
    std::string ID;
    int stereo;
    int satellite;
    std::string satellite_location;
    std::string request_location;
    std::string time;
    double area;
    int strips;
    double duration;
    double distance;
    double angle;
    double sun_elevation;
    double cloud_cover_estimate;
    int priority;
    int priority_mod;
    int customer_type_mod;
    double price;
    double waiting_time;
    double uncertainty;
    double cloud_cover_real;
    double score_scenario;
    std::string score_method;
    double score_alpha;

    // Pre-computed fields (populated in ReadData, avoid repeated parsing)
    long long epoch_seconds = 0;
    Vec3 sat_xyz = {0, 0, 0};
    Vec3 req_xyz = {0, 0, 0};
};

struct InstanceData
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
};

// -----------------------------
// Utilitários
// -----------------------------
static inline std::string trim(std::string s) {
    auto not_space = [](unsigned char ch){ return !std::isspace(ch); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

static std::vector<std::string> split_csv_line(const std::string& line)
{
    std::vector<std::string> fields;
    std::string cur;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i)
    {
        char c = line[i];

        if (c == '"')
        {
            in_quotes = !in_quotes; // alterna estado
            continue;               // não guarda aspas
        }

        if (c == ',' && !in_quotes)
        {
            fields.push_back(trim(cur));
            cur.clear();
        }
        else
        {
            cur.push_back(c);
        }
    }

    fields.push_back(trim(cur));
    return fields;
}

static inline bool starts_with(const std::string& s, const std::string& prefix) {
    return s.size() >= prefix.size() && s.compare(0, prefix.size(), prefix) == 0;
}

static std::vector<int> parse_int_list_in_brackets(const std::string& line) {
    // Espera algo como: "... [1, 2, 3]" ou "... [38755]"
    auto l = line.find('[');
    auto r = line.find(']');
    if (l == std::string::npos || r == std::string::npos || r <= l) {
        throw std::runtime_error("Lista entre colchetes não encontrada: " + line);
    }

    std::string inside = line.substr(l + 1, r - l - 1);
    inside = trim(inside);

    std::vector<int> out;
    if (inside.empty()) return out;

    std::stringstream ss(inside);
    std::string token;
    while (std::getline(ss, token, ',')) {
        token = trim(token);
        if (!token.empty()) {
            out.push_back(std::stoi(token));
        }
    }
    return out;
}

static DateTime parse_datetime_list(const std::string& line) {
    // Espera: "Início: [2025, 11, 21, 9, 40]"
    auto vals = parse_int_list_in_brackets(line);
    if (vals.size() != 5) {
        throw std::runtime_error("Data/hora deve ter 5 inteiros [Y,M,D,H,m]: " + line);
    }
    DateTime dt;
    dt.year = vals[0];
    dt.month = vals[1];
    dt.day = vals[2];
    dt.hour = vals[3];
    dt.minute = vals[4];

    if (!dt.valid()) {
        throw std::runtime_error("Data/hora inválida em: " + line);
    }
    return dt;
}

static int parse_horizon_hours(const std::string& line) {
    // Espera: "Horizonte: 8 horas" (também tolera "hora" singular)
    auto pos = line.find(':');
    if (pos == std::string::npos) {
        throw std::runtime_error("Formato inválido para horizonte: " + line);
    }
    std::string after = trim(line.substr(pos + 1));
    // Pega o primeiro inteiro que aparecer
    std::stringstream ss(after);
    int h = 0;
    ss >> h;
    if (!ss || h <= 0) {
        throw std::runtime_error("Horizonte (horas) inválido em: " + line);
    }
    return h;
}

// -----------------------------
// Leitura do arquivo info
// -----------------------------
InstanceData read_eos_instance(const std::string& info_path)
{
    InstanceData data;

    // =========================
    // 1) Ler info.txt
    // =========================
    std::ifstream fin(info_path);
    if (!fin.is_open())
        throw std::runtime_error("Não foi possível abrir o info.txt: " + info_path);

    std::string line;

    while (std::getline(fin, line))
    {
        line = trim(line);

        if (starts_with(line, "Início:"))
        {
            DateTime dt = parse_datetime_list(line);
            data.start_year = dt.year;
            data.start_month = dt.month;
            data.start_day = dt.day;
            data.start_hour = dt.hour;
            data.start_minute = dt.minute;
        }
        else if (starts_with(line, "Horizonte:"))
        {
            data.horizon_hours = parse_horizon_hours(line);
        }
        else if (starts_with(line, "NORADs usados:"))
        {
            data.norads = parse_int_list_in_brackets(line);
        }
    }

    fin.close();

    // =========================
    // 2) Descobrir caminho do CSV
    // =========================
    std::string csv_path = info_path;
    size_t pos = csv_path.find("_info.txt");
    if (pos != std::string::npos)
        csv_path.replace(pos, 9, "_pf_df.csv");
    else
        throw std::runtime_error("Nome do arquivo info.txt não contém '_info.txt': " + info_path);

    // =========================
    // 3) Ler CSV
    // =========================
    std::ifstream csv(csv_path);
    if (!csv.is_open())
        throw std::runtime_error("Não foi possível abrir o pf_df.csv: " + csv_path);

    // pular header
    std::getline(csv, line);

    data.acquisitions.clear();

    while (std::getline(csv, line))
    {
        auto cols = split_csv_line(line);

        // Esperamos 24 colunas, conforme seu header
        if (cols.size() < 24)
        {
            std::cout << "Linha CSV com colunas insuficientes (" << cols.size() << "):\n"
                      << line << "\n";
            continue;
        }

        Acquisition acq;
        acq.index = std::stoi(cols[0]);
        acq.ID = cols[1];
        acq.stereo = std::stoi(cols[2]);
        acq.satellite = std::stoi(cols[3]);
        acq.satellite_location = cols[4];
        acq.request_location = cols[5];
        acq.time = cols[6];
        acq.area = std::stod(cols[7]);
        acq.strips = std::stoi(cols[8]);
        acq.duration = std::stod(cols[9]);
        acq.distance = std::stod(cols[10]);
        acq.angle = std::stod(cols[11]);
        acq.sun_elevation = std::stod(cols[12]);
        acq.cloud_cover_estimate = std::stod(cols[13]);
        acq.priority = std::stoi(cols[14]);
        acq.priority_mod = std::stoi(cols[15]);
        acq.customer_type_mod = std::stoi(cols[16]);
        acq.price = std::stod(cols[17]);
        acq.waiting_time = std::stod(cols[18]);
        acq.uncertainty = std::stod(cols[19]);
        acq.cloud_cover_real = std::stod(cols[20]);
        acq.score_scenario = std::stod(cols[21]);
        acq.score_method = cols[22];
        acq.score_alpha = std::stod(cols[23]);

        data.acquisitions.push_back(acq);
    }

    csv.close();

    data.n = static_cast<int>(data.acquisitions.size());

    return data;
}
