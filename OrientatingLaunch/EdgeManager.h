#pragma once
// EdgeManager.h
#ifndef EDGEMANAGER_HPP
#define EDGEMANAGER_HPP

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <memory>

    // Ç°ÖÃÉùÃ÷
class EdgeManager;

class DirectedEdge {
public:
    std::string from;
    std::string to;
    double azimuth;
    DirectedEdge* reverse_edge;
    bool isFixed;

    DirectedEdge(const std::string& f, const std::string& t, double az, bool isFixed = false);
    void setAzimuth(double new_azimuth);
    bool isAzimuthFixed() const { return isFixed; }

private:
    void setAzimuthInternal(double new_azimuth, bool propagate);
    friend class EdgeManager;
};

class EdgeManager {
private:
    struct PairHash {
        template <class T1, class T2>
        std::size_t operator () (const std::pair<T1, T2>& p) const {
            auto h1 = std::hash<T1>{}(p.first);
            auto h2 = std::hash<T2>{}(p.second);
            return h1 ^ (h2 << 1);
        }
    };

    std::unordered_map<std::pair<std::string, std::string>, DirectedEdge*, PairHash> edge_map;

public:
    using ConstIterator = decltype(edge_map)::const_iterator;

    EdgeManager();
    ~EdgeManager();

    DirectedEdge* addEdge(const std::string& from, const std::string& to,
        double azimuth, bool isFixed = false);
    DirectedEdge* getEdge(const std::string& from, const std::string& to) const;
    void updateAzimuth(const std::string& from, const std::string& to, double delta);
    void printEdges() const;

    ConstIterator begin() const { return edge_map.begin(); }
    ConstIterator end() const { return edge_map.end(); }
};





#endif // EDGEMANAGER_HPP