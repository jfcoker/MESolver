#include "pch.h"
#include "site.h"
#include "consts.h"

site::site(double x, double y, double z, double E) 
{
    pos.X = x;
    pos.Y = y;
    pos.Z = z;
    energy = E;
}

void site::addNeighbour(site* pSite, double J)
{
    // Unless already present, add passed site as neighbour to 'this'
    if (!this->hasNeighbour(pSite))
    {
        neighbour* n = new site::neighbour;
        n->_pSite = pSite;
        n->_J = J;
        neighbours.push_back(n);
    }

    // Unless already present, add 'this' as neighbour to passed site
    if (!pSite->hasNeighbour(this))
    {
        neighbour* th = new site::neighbour;
        th->_pSite = this;
        th->_J = J;
        pSite->neighbours.push_back(th);
    }

}

site::neighbour* site::hasNeighbour(site* pSite)
{
    std::vector<neighbour*>::iterator it = neighbours.begin();
    for (int i = 0; it != neighbours.end(); i++, it++)
        if ((*it)->_pSite == pSite)
            return *it;

    return NULL;
}

double site::deltaE(site* pSite, double fieldZ)
{
    double deltaE = (pSite->energy - this->energy) + (pSite->pos.Z - this->pos.Z) * fieldZ; // deltaE = (ep_f - ep_i) + e * [M_f.Z - M_i.Z] * F_z
}

double site::Rate(site* pSite, double fieldZ, double kBT, double reorg)
{
    neighbour* pN = this->hasNeighbour(pSite);

    if (!pN)
        // Sites aren't interacting, transfer rate will be zero.
        return 0.0;
    else
    {
        double J2 = std::pow(pN->_J,2); // |J_if|^2
        double deltaE = this->deltaE(pSite, fieldZ);
        return ((2 * pi) / hbar) * J2 * std::pow(4 * pi * reorg * kBT, -0.5) * std::exp(-1 * std::pow(deltaE + reorg, 2) / (4 * reorg * kBT));
    }
}

site::~site()
{
    std::vector<neighbour*>::iterator it = neighbours.begin();
    for (int i = 0; it != neighbours.end(); i++, it++)
        delete (*it);
    
    neighbours.clear();
}


std::ostream& operator<<(std::ostream& os, const site& st)
{
    os << "pos = (" << st.pos.X << "," << st.pos.Y << "," << st.pos.Z << "); E = " << st.energy << "; # neighbors = " << st.neighbours.size();
    return os;
}