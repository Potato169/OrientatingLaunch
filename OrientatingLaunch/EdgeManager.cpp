// EdgeManager.cpp
#include "EdgeManager.h"



// DirectedEdge 成员函数实现
DirectedEdge::DirectedEdge(const std::string& f, const std::string& t,
    double az, bool isFixed)
    : from(f), to(t), azimuth(az), reverse_edge(nullptr), isFixed(isFixed) {
    /*setAzimuthInternal(az, false);*/
    if (az < 0.0 || az >= 360.0) {
        throw std::invalid_argument("Azimuth must be in [0, 360) degrees");
    }
    if (reverse_edge) {
        reverse_edge->setAzimuthInternal(fmod(az + 180.0, 360.0), false);
        /*reverse_edge->isFixed = isFixed;*/
    }
}

void DirectedEdge::setAzimuth(double new_azimuth) {
    setAzimuthInternal(new_azimuth, true);
}



void DirectedEdge::setAzimuthInternal(double new_azimuth, bool propagate) {
    if (isFixed) {
        throw std::invalid_argument("Azimuth of fixed edge cannot be modified.");
    }
    if (new_azimuth < 0.0 || new_azimuth >= 360.0) {
        throw std::invalid_argument("Azimuth must be in [0, 360) degrees");
    }
    azimuth = new_azimuth;
    if (propagate && reverse_edge) {
        reverse_edge->setAzimuthInternal(fmod(new_azimuth + 180.0, 360.0), false);
    }
}

// EdgeManager 成员函数实现
EdgeManager::EdgeManager() = default;

EdgeManager::~EdgeManager() {
    std::unordered_set<DirectedEdge*> deleted_edges;
    for (auto& pair : edge_map) {
        if (deleted_edges.find(pair.second) == deleted_edges.end()) {
            delete pair.second;
            deleted_edges.insert(pair.second->reverse_edge);
        }
    }
}

DirectedEdge* EdgeManager::addEdge(const std::string& from, const std::string& to,
    double azimuth, bool isFixed) {
    auto key = std::make_pair(from, to);
    auto reverse_key = std::make_pair(to, from);

    bool new_exists = edge_map.find(key) != edge_map.end();
    bool reverse_exists = edge_map.find(reverse_key) != edge_map.end();

    try {
        if (new_exists || reverse_exists) {
            if (new_exists) {
                DirectedEdge* existing = edge_map.at(key);
                std::cout << "Err: Edge " << from << "→" << to
                    << " already exists with azimuth "
                    << std::fixed << std::setprecision(8) << existing->azimuth << "°\n";
                return existing;
            }
            else {
                DirectedEdge* reverse_edge = edge_map.at(reverse_key);
                std::cout << "Err: Edge " << from << "→" << to
                    << " already exists via reverse edge "
                    << reverse_edge->from << "→" << reverse_edge->to
                    << " with azimuth " << std::fixed << std::setprecision(8)
                    << reverse_edge->azimuth << "°\n";
                return reverse_edge->reverse_edge;
            }
        }

        double reverse_az = fmod(azimuth + 180.0, 360.0);
        if (reverse_az < 0) reverse_az += 360.0;

        DirectedEdge* new_edge = new DirectedEdge(from, to, azimuth, isFixed);
        DirectedEdge* reverse_edge = new DirectedEdge(to, from, reverse_az, isFixed);
        new_edge->reverse_edge = reverse_edge;
        reverse_edge->reverse_edge = new_edge;

        edge_map.emplace(key, new_edge);
        edge_map.emplace(reverse_key, reverse_edge);

        // 更新计数器
        if (isFixed) {
            fixed_edge_count++;
        }

        std::cout << " Created edge " << from << "→" << to
            << " with azimuth " << std::fixed << std::setprecision(8) << azimuth << "°\n"
            << "and reverse edge " << to << "→" << from
            << " with azimuth " << std::fixed << std::setprecision(8) << reverse_az << "°\n";

        return new_edge;

    }
    catch (const std::exception& e) {
        std::cerr << " Failed to create edge " << from << "→" << to
            << ": " << e.what() << "\n";
        return nullptr;
    }
}

DirectedEdge* EdgeManager::getEdge(const std::string& from, const std::string& to) const {
    auto it = edge_map.find(std::make_pair(from, to));
    return it != edge_map.end() ? it->second : nullptr;
}

void EdgeManager::updateAzimuth(const std::string& from, const std::string& to, double delta) {
    if (DirectedEdge* edge = getEdge(from, to)) {
        edge->setAzimuth(edge->azimuth + delta);
    }
    else if (DirectedEdge* revEdge = getEdge(to, from)) {
        revEdge->setAzimuth(revEdge->azimuth - delta);
    }
}

void EdgeManager::printEdges() const {
    std::unordered_set<const DirectedEdge*> visited;

    std::cout << "----- Edge List -----\n";
    for (const auto& pair : edge_map) {
        const DirectedEdge* edge = pair.second;

        if (visited.find(edge) != visited.end()) continue;

        std::ostringstream oss1, oss2;
        oss1 << std::left << std::setw(6) << edge->from << "→" << std::setw(6) << edge->to
            << " (Azimuth: " << std::fixed << std::setprecision(8) << edge->azimuth << "°)";

        oss2 << std::left << std::setw(6) << edge->reverse_edge->from << "→"
            << std::setw(6) << edge->reverse_edge->to
            << " (Azimuth: " << std::fixed << std::setprecision(8)
            << edge->reverse_edge->azimuth << "°)";

        std::cout << "Edge " << oss1.str() << " | Reverse Edge " << oss2.str() << "\n";

        visited.insert(edge);
        visited.insert(edge->reverse_edge);
    }
    std::cout << "---------------------\n";
}

