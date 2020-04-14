#include "pch.h"
#include "site.h"

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

double site::Rate(site* pSite)
{
    if (!this->hasNeighbour(pSite))
        return 0.0;
    else
        return 1.0;

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