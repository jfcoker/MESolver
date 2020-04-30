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
    return (pSite->energy - this->energy) + (pSite->pos.Z - this->pos.Z) * fieldZ; // deltaE = (ep_f - ep_i) + [M_f.Z - M_i.Z] * F_z
}

double site::Rate(site* pSite, double fieldZ, double kBT, double reorg)
{
    neighbour* pN = this->hasNeighbour(pSite);

    if (!pN)
        // Sites aren't interacting, transfer rate will be zero.
        return 0.0;
    else if (pN->_rate < 0.0) // If rate is -1, then this is the initial call and rate needs to be calculated. Using '<' to avoid equality comparison on floating point values.
    {
        double J2 = std::pow(pN->_J, 2); // |J_if|^2
        double deltaE = this->deltaE(pSite, fieldZ);
        pN->_rate = ((2 * pi) / hbar) * J2 * std::pow(4 * pi * reorg * kBT, -0.5) * std::exp(-1 * std::pow(deltaE + reorg, 2) / (4 * reorg * kBT));
        return pN->_rate;
    }
    else // Rate has already been calculated. Just grab it.
        return pN->_rate;
}

double site::PrecondFactor(double fieldZ, double kBT, double reorg, double E0, site::PrecondForm form, bool apply)
{
    if (apply) 
    {
        // Could use interface / implementation instead of enum / switch
        switch (form)
        {
        case PrecondForm::boltzmann: 
            return std::exp((E0 - (this->energy + this->pos.Z * fieldZ)) / kBT);

        case PrecondForm::boltzmannSquared: 
            return std::pow(std::exp((E0 - (this->energy + this->pos.Z * fieldZ)) / kBT), 2.0);

            //ALTERNATIVE PRCONDITIONING FACTOR FROM: Yu PRB 2001 (long)
        case PrecondForm::rateSum:
        {
            double sum = 0;
            std::vector<neighbour*>::iterator it = neighbours.begin();
            for (int i = 0; it != neighbours.end(); i++, it++)
                sum += (*it)->_pSite->Rate(this,fieldZ,kBT,reorg);

            return sum;
        }
         default:
             throw std::logic_error("Form of preconditioning factor not implemented.");
             return -1.0;

        }
    }
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