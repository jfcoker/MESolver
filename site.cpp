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
        neighbour* n = new neighbour;
        n->_pSite = pSite;
        n->_J = J;
        neighbours.push_back(n);
    }

    // Unless already present, add 'this' as neighbour to passed site
    if (!pSite->hasNeighbour(this))
    {
        neighbour* th = new neighbour;
        th->_pSite = this;
        th->_J = J;
        pSite->neighbours.push_back(th);
    }

}

neighbour* site::hasNeighbour(site* pSite)
{
    std::vector<neighbour*>::iterator it = neighbours.begin();
    for (int i = 0; it != neighbours.end(); i++, it++)
        if ((*it)->_pSite == pSite)
            return *it;

    return NULL;
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