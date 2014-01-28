#include "interaction.hpp"

#include <vector>
//#include <set>
//#include <map>
//#include <cassert>
//#include <climits>

#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

class constraint {
    public:
        virtual void apply(Box &box) = 0;
        virtual int ndof() = 0;
        virtual ~constraint(){};
};

class coordConstraint : public constraint {
    private:
        atom* a;
        bool fixed[3];
        Vec loc;
    public:
        coordConstraint(atom* atm, bool fixx, bool fixy, bool fixz, Vec loc) :
            a(atm), loc(loc) {fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        coordConstraint(atom* atm, bool fixx, bool fixy, bool fixz) :
            a(atm), loc(a->x) {fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        coordConstraint(atom* atm) :
            a(atm), loc(a->x) {fixed[0] = fixed[1] = fixed[2] = true;};
        int ndof(){return (int)fixed[0] + (int)fixed[1] + (int)fixed[2];};
        void apply(Box &box){
            for(uint i=0; i<3; i++){
                if(not fixed[i]) continue;
                a->f[i] = 0;
                a->v[i] = 0;
                a->x[i] = loc[i];
            }
        }
};

class coordCOMConstraint : public constraint {
    private:
        atomgroup* a;
        bool fixed[3];
        Vec loc;
    public:
        coordCOMConstraint(atomgroup* atm, bool fixx, bool fixy, bool fixz, Vec loc) :
            a(atm), loc(loc) {fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        coordCOMConstraint(atomgroup* atm, bool fixx, bool fixy, bool fixz) :
            a(atm), loc(a->com()) {fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        coordCOMConstraint(atomgroup* atm) :
            a(atm), loc(a->com()) {fixed[0] = fixed[1] = fixed[2] = true;};
        int ndof(){return (int)fixed[0] + (int)fixed[1] + (int)fixed[2];};
        void apply(Box &box){
            Vec com = a->com() - loc;
            Vec comv = a->comv();
            Vec totf = Vec();
            for(uint i=0; i< a->size(); i++){
                totf += (*a)[i].f;
            }
            Vec tota = totf / a->mass();
            
            for(uint i=0; i< a->size(); i++){
                atom &atm = (*a)[i];
                Vec df = (tota * (atm.m));
                for(uint j=0; j<3; j++){
                    if(not fixed[j]) continue;
                    atm.f[j] -= df[j];
                    atm.v[j] -= comv[j];
                    atm.x[j] -= com[j];
                }
            }
        }
};

class relativeConstraint : public constraint {
    private:
        atom *a1, *a2;
        bool fixed[3];
        Vec loc;
    public:
        relativeConstraint(atom* atm1, atom* atm2, bool fixx, bool fixy, bool fixz, Vec loc) :
            a1(atm1), a2(atm2), loc(loc) {
                fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        relativeConstraint(atom* atm1, atom* atm2, bool fixx, bool fixy, bool fixz) :
            a1(atm1), a2(atm2), loc(a2->x - a1->x) {
                fixed[0] = fixx; fixed[1] = fixy; fixed[2] = fixz;};
        relativeConstraint(atom* atm1, atom* atm2) :
            a1(atm1), a2(atm2), loc(a2->x - a1->x) {
                fixed[0] = fixed[1] = fixed[2] = true;};
        int ndof(){return (int)fixed[0] + (int)fixed[1] + (int)fixed[2];};
        void apply(Box &box){
            flt mratio1 = a1->m / (a1->m + a2->m);
            flt mratio2 = a2->m / (a1->m + a2->m);
            Vec totf = a2->f + a1->f;
            Vec dv = a2->v - a1->v;
            Vec dx = a2->x - a1->x;
            for(uint i=0; i<NDIM; i++){
                if(not fixed[i]) continue;
                a1->f[i] = totf[i]*mratio1;
                a2->f[i] = totf[i]*mratio2;
                a1->v[i] += dv[i]/2;
                a2->v[i] -= dv[i]/2;
                a1->x[i] += dx[i]/2;
                a2->x[i] -= dx[i]/2;
                assert(abs(a2->v[i] - a1->v[i]) < 1e-5);
                assert(abs(a2->x[i] - a1->x[i]) < 1e-5);
            }
        }
};

class distConstraint : public constraint {
    private:
        atom *a1, *a2;
        flt dist;
    public:
        distConstraint(atom* atm1, atom* atm2, flt dist) :
            a1(atm1), a2(atm2), dist(dist) {};
        int ndof(){return 1;};
        void apply(Box &box){
            flt M = (a1->m + a2->m);
            flt mratio1 = a1->m / M;
            flt mratio2 = a2->m / M;
            
            Vec dx = a2->x - a1->x;
            flt dxmag = dx.mag();
            Vec dxnorm = dx / dxmag;
            
            a1->x += dx * ((1 - dist/dxmag)*mratio2);
            a2->x -= dx * ((1 - dist/dxmag)*mratio1);
            //~ dx = a2->x - a1->x;
            //~ dxmag = dx.mag();
            //~ dxnorm = dx / dxmag; dxnorm should still be the same
            
            Vec baddv = dxnorm * ((a2->v - a1->v).dot(dxnorm)/2);
            a1->v += baddv;
            a2->v -= baddv;
            
            // newv2 • u = v2 - |baddv|
            // newv1 • u = v1 + |baddv|
            // (newv2 - newv1) • u = (v2 - v1 - (2*baddv)) • u
            //                     = ((v2 - v1)•u - (2*baddv)•u)
            //                     = (|baddv|*2 - |baddv|*2) = 0
            assert((a2->v - a1->v).dot(dxnorm) < 1e-8);
            
            // TODO: Fix mass ratio stuff
            Vec baddf = dxnorm * ((a2->f - a1->f).dot(dxnorm)/2);
            a1->f += baddf;
            a2->f -= baddf;
            assert((a2->f - a1->f).dot(dxnorm) < 1e-8);
        }
};


class linearConstraint : public constraint {
    private:
        atomgroup& atms;
        flt dist;
        flt lincom, I, M;
    public:
        linearConstraint(atomgroup& atms, flt dist) :
            atms(atms), dist(dist), lincom(0), I(0), M(0) {
            for(uint i = 0; i < atms.size(); i++){
                M += atms[i].m;
                lincom += (dist*i)*atms[i].m;
            }
            lincom /= M;
            
            for(uint i = 0; i < atms.size(); i++){
                flt dx = (dist*i - lincom);
                I += atms[i].m * dx * dx;
            }
        };
        int ndof(){return atms.size()-1;};
        
        void apply(Box &box){
            Vec com = atms.com();
            Vec comv = atms.comv();
            Vec comf = Vec();
            
            uint sz = atms.size();
            Vec lvec = Vec();
            #ifdef VEC3D
            Vec L = Vec();
            Vec omega  = Vec();
            Vec tau = Vec();
            Vec alpha = Vec();
            #else
            flt L = 0;
            flt omega = 0;
            flt tau = 0;
            flt alpha = 0;
            #endif
            
            for(uint i = 0; i < sz; i++){
                flt chaindist = i * dist - lincom;
                Vec dx = atms[i].x - com;
                comf += atms[i].f;
                lvec += dx.norm() * chaindist;
                L += dx.cross(atms[i].v) * atms[i].m;
                tau += dx.cross(atms[i].f);
            }
            
            lvec.normalize();
            omega = L / I;
            alpha = tau / I;
            
            for(uint i = 0; i < sz; i++){
                flt chaindist = i * dist - lincom;
                Vec dx = lvec*chaindist;
                atms[i].x = com + dx;
                atms[i].v = comv + dx.cross(omega);
                atms[i].f = comf + (dx.cross(alpha)*atms[i].m);
            }
        }
};

class NPHGaussianConstraint : public constraint {
    private:
        sptr<OriginBox> box;
        flt ddV, dV; // that's dV²/dt², dV/dt
        vector<atomgroup*> groups;
    public:
        NPHGaussianConstraint(sptr<OriginBox> box, vector<atomgroup*> groups) : 
                box(box), ddV(0), dV(0), groups(groups){};
        int ndof(){return 0;};
        void apply(Box &box2){
            //~ flt V = box->V();
            assert((Box*) box.get() == &box2);
            //~ vector<atomgroup*>::iterator git;
            //~ for(git = groups.begin(); git<groups.end(); git++){
                //~ atomgroup &m = **git;
                //~ for(uint i=0; i<m.size(); i++){
                    //~ m[i].v += 
                //~ }
            //~ }
        };
};
#endif


