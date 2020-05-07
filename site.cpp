#include "pch.h"
#include "site.h"
#include "transporter.h"
#include "consts.h"
#include "IO.h"

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

site::neighbour* site::hasNeighbour(site* pSite)
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

std::vector<site> CreateSites(char* XYZfile, char* EDGEfile)
{
    std::vector<site> sites;

    // Use contents of .xyz file to populate site vector
    // Columns should be formatted as: x (float), y (float), z (float), site type (char, unused), site energy (float)
    std::ifstream in;
    open(XYZfile, in);
    std::string line;

    int linenum = 0;
    while (std::getline(in, line))
    {
        linenum++;
        std::stringstream linestream(line);
        double x, y, z, E;
        char s;
        linestream >> x >> y >> z >> s >> E;

        if (linestream.fail())
        {
            std::cout << "***ERROR***: Unexpected formatting of .xyz file at line " << linenum << ". Expected 5 columns with value types; float, float, float, char, float.\n";
            in.close();
            exit(-1);
        }

        sites.push_back(site(x, y, z, E));
    }

    in.close();

    // Use contents of .edge file to set interacting neighbours
    // Columns should be formatted as: site 1 (int), site 2 (int), J (float)
    open(EDGEfile, in);

    linenum = 0;
    while (std::getline(in, line))
    {
        linenum++;
        std::stringstream linestream(line);
        int s1, s2;
        double J;
        linestream >> s1 >> s2 >> J;

        if (linestream.fail())
        {
            std::cout << "***ERROR***: Unexpected formatting of .edge file at line " << linenum << ". Expected 3 columns with value types; int, int, float.\n";
            in.close();
            exit(-1);
        }

        if (s1 >= sites.size() || s2 >= sites.size())
        {
            std::cout << "***ERROR***: Site index out of range in .edge file at line " << linenum << "\n";
            in.close();
            exit(-1);
        }

        sites[s1].addNeighbour(&sites[s2], J);

    }

    in.close();

    return sites;

}