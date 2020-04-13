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
    bool isPassedPresent = false;
    std::vector<neighbour*>::iterator it1 = neighbours.begin();
    for (int i = 0; it1 != neighbours.end(); i++, it1++) 
        if ((*it1)->_pSite == pSite)
            isPassedPresent = true;
    
    if (!isPassedPresent)
    {
        neighbour* n = new site::neighbour;
        n->_pSite = pSite;
        n->_J = J;
        neighbours.push_back(n);
    }


    // Unless already present, add 'this' as neighbour to passed site
    bool isThisPresent = false;
    std::vector<neighbour*>::iterator it2 = pSite->neighbours.begin();
    for (int i = 0; it2 != pSite->neighbours.end(); i++, it2++) 
        if ((*it2)->_pSite == this)
            isThisPresent = true;
    
    if (!isThisPresent) 
    {
        neighbour* th = new site::neighbour;
        th->_pSite = this;
        th->_J = J;
        pSite->neighbours.push_back(th);
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