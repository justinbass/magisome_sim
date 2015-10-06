#include <cmath>
#include <xpre.h>
using namespace std;

#define timestep 1

const double pi = 3.14152653589793238;
const double hbar = 6.6260695729e-34;
const double e = 1.60217656535e-19;
const double c = 2.9979e8;
const double amu = 1.66053892172e-27;

class Photon;
class Positron;
class Neutrino;
class Nucleon;

class Particle {
public:
	//Particle neighbors
	Particle* right;
	Particle* left;

	//Particle properties:
	double E; //Energy
	double M; //Mass
	double q; //Charge
	double x; //Position
    const char* name;

	//Mass (amu)
	double get_u() {
		return M / amu;
	}

	void set_u(double uu) {
		M = uu * amu;
	}

	//Energy (eV)
	double get_eV() {
		return this->E / e;
	}

	double set_eV(double ee) {
		this->E  = ee / e;
	}

	//Velocity (m/s)
	double get_v() {
		struct xpr xE = dbltox(E);
		struct xpr xM = dbltox(M);
		struct xpr xCC = dbltox(c*c);

		return c*sqrt(xtodbl(xadd(xOne, xdiv(xOne, xpwr(xadd(xOne, xdiv(xE, xmul(xM, xCC)), 0), 2)), 1)));
		//return c*sqrt(1.0 - (1.0 / ((1.0 + (E / (M*c*c)))*(1.0 + (E / (M*c*c))))));
	}

    //Remove 
    double remove() {
        if (right and right->left) {
            printf("RIGHT %s\n", right->name);
            printf("RIGHT LEFT %s\n", right->left->name);
            right->left = left;
        }

        if (left)
            left->right = right;
    }

    void add_to_right(Particle* p) {
        p->right = this->right
        this->right = p
    }

    void add_to_left(Particle* p) {
        p->left = this->left
        this->left = p
    }

	//Particle interactions, abstract methods to be defined in the subclass
	virtual void interact(Particle* other) = 0;
	virtual void interact(Photon* other) = 0;
	virtual void interact(Positron* other) = 0;
	virtual void interact(Neutrino* other) = 0;
	virtual void interact(Nucleon* other) = 0;
	
};

class Particles {
public:
    Particle* middle;

    Particle* get_leftmost() {
        //Get leftmost particle
        Particle* pleft = middle;
        while(pleft->left)
            pleft = pleft->left;

        return pleft;
    }


    int check_validity() {
        Particle *p1 = get_leftmost();
        Particle *p2 = p1->right;

        //Make sure particles to the right have increasing position
        while(p2) {
            if(p1->x >= p2->x)
                return false;
            p1 = p1->right;
            p2 = p2->right;
        }

        return true;

    }

    void print() {
        printf("Particles:\n");

        Particle* pleft = get_leftmost();

        //Print all particles
        while(pleft) {
            if(pleft == middle)
                printf("*");

            printf("%s at %lf nm\n", pleft->name, pleft->x*1e9);
            pleft = pleft->right;
        }

        return;
    }

    void interact() {
        printf("\nInteractions:\n");
        int seconds = 1;
        int iterations_per_second = 1;
        for (int i = 0; i < iterations_per_second; i++) {
            Particle* current = get_leftmost();
            while (current->right != NULL) {
                current->interact(current->right);
                current = current->right;
            }
        }
    }

    Particles(Particle *mmiddle) {
        middle = mmiddle;
    }
};


class Products {
public:
	//The left and right-most products of the Products. Note that the two
	//do not have to be neighbors: several products can exist between
	Particle* product_left;
	Particle* product_right;
};

class Photon : public Particle {
public:
	Photon(double EE, double xx) {
		M = 0.0;
		q = 0.0;
		E = EE;
		x = xx;
        name = "Photon";
	}

	virtual void interact(Particle* other);
	virtual void interact(Photon* other);
	virtual void interact(Positron* other);
	virtual void interact(Neutrino* other);
	virtual void interact(Nucleon* other);
};

class Positron : public Particle {
public:
	Positron(double EE, double xx) {
		M = 9.1093829140e-31;
		q = 1.0;
		E = EE;
		x = xx;
        name = "Positron";
	}

	virtual void interact(Particle* other);
	virtual void interact(Photon* other);
	virtual void interact(Positron* other);
	virtual void interact(Neutrino* other);
	virtual void interact(Nucleon* other);
};

class Neutrino : public Particle {
public:
	Neutrino(double EE, double xx) {
		set_u(0.2);
		q = 0.0;
		E = EE;
		x = xx;
        name = "Neutrino";
	}

	virtual void interact(Particle* other);
	virtual void interact(Photon* other);
	virtual void interact(Positron* other);
	virtual void interact(Neutrino* other);
	virtual void interact(Nucleon* other);
};

//Nucleons
class Proton;
class Deuteron;
class Triton;
class Helion3;
class Helion4;

class Nucleon : public Particle {
public:
	int Z; //Number of protons
	int N; //Number of nucleons
	double r; //Radius (m)
	
	//Gamow Probability - the probability of barrier penetration for two Nucleons
	double gamow_probability(Nucleon* b){
		double reduced_mass = this->get_u() * b->get_u() / (this->get_u() + b->get_u());
		return exp(-pi * e * e * this->Z * b->Z * sqrt(2 * reduced_mass / this->E) / hbar);
	}

	//Cross-section in m^2
	double cross_section(Nucleon *b, double SE) {
		return SE * this->gamow_probability(b) / this->E;
	}

	//Total probability of fusion
	double fusion_probability(Nucleon *b, double SE) {
		return pi * this->r * this->r / this->cross_section(b, SE);
	}

	virtual void interact(Particle* other);
	virtual void interact(Photon* other);
	virtual void interact(Positron* other);
	virtual void interact(Neutrino* other);
	virtual void interact(Nucleon* other);

	virtual Products* products(Nucleon* other) = 0;
	virtual Products* products(Proton* other) = 0;
	virtual Products* products(Deuteron* other) = 0;
	virtual Products* products(Triton* other) = 0;
	virtual Products* products(Helion3* other) = 0;
	virtual Products* products(Helion4* other) = 0;
};

class Proton : public Nucleon {
public:
	Proton(double EE, double xx) {
		Z = q = 1;
		N = 0;
		set_u(1.0072765);
		r = 0.878e-15;
		E = EE;
		x = xx;
        name = "Proton";
	}

	//NOTE: Products returned here, and random generation done here as well
	virtual Products* products(Nucleon* other) { return other->products(this); }
	virtual Products* products(Proton* other) {
		double SE = 6.41e-66;
		Products* ProtonProtonProducts = new Products();
		fusion_probability(other, SE);
		
		return ProtonProtonProducts;
	}
	virtual Products* products(Deuteron* other) { double SE = 4.0e-48; return NULL; }
	virtual Products* products(Triton* other) { double SE = r*r; return NULL; }
	virtual Products* products(Helion3* other) { double SE = 3.7e-64; return NULL; }
	virtual Products* products(Helion4* other) { return NULL; }

};

class Deuteron : public Nucleon {
public:
	Deuteron(double EE, double xx) {
		Z = q = 1;
		N = 1;
		set_u(2.0135532);
		r = 2.14e-15;
		E = EE;
		x = xx;
        name = "Deuteron";
	}

	virtual Products* products(Nucleon* other) { return other->products(this); }
	virtual Products* products(Proton* other) { double SE = 4.0e-48; return NULL; }
	virtual Products* products(Deuteron* other) { double SE = r*r; return NULL; }
	virtual Products* products(Triton* other) { double SE = r*r; return NULL; }
	virtual Products* products(Helion3* other) { double SE = 9.5e-41; return NULL; }
	virtual Products* products(Helion4* other) { return NULL; }

};

class Triton : public Nucleon {
public:
	Triton(double EE, double xx) {
		Z = q = 1;
		N = 2;
		set_u(3.015501);
		r = 1.76e-15;
		E = EE;
		x = xx;
        name = "Triton";
	}

	virtual Products* products(Nucleon* other) { return other->products(this); }
	virtual Products* products(Proton* other) { double SE = r*r; return NULL; }
	virtual Products* products(Deuteron* other) { double SE = r*r; return NULL; }
	virtual Products* products(Triton* other) { double SE = r*r; return NULL; }
	virtual Products* products(Helion3* other) { double SE = r*r; return NULL; }
	virtual Products* products(Helion4* other) { return NULL; }

};

class Helion3 : public Nucleon {
public:
	Helion3(double EE, double xx) {
		Z = q = 2;
		N = 1;
		set_u(3.015481);
		r = 1.97e-15;
		E = EE;
		x = xx;
        name = "Helion3";
	}

	virtual Products* products(Nucleon* other) { other->products(this); }
	virtual Products* products(Proton* other) { double SE = 3.7e-64; return NULL; }
	virtual Products* products(Deuteron* other) { double SE = 9.5e-41; return NULL; }
	virtual Products* products(Triton* other) { double SE = r*r; return NULL; }
	virtual Products* products(Helion3* other) { double SE = 8.7e-41; return NULL; }
	virtual Products* products(Helion4* other) { return NULL; }

};

class Helion4 : public Nucleon {
public:
	Helion4(double EE, double xx) {
		Z = q = 2;
		N = 2;
		set_u(4.00150618);
		r = 1.68e-15;
		E = EE;
		x = xx;
        name = "Helion4";
	}

	virtual Products* products(Nucleon* other) { other->products(this); }
	virtual Products* products(Proton* other) { double SE = 3.7e-64; return NULL; }
	virtual Products* products(Deuteron* other) { double SE = 9.5e-41; return NULL; }
	virtual Products* products(Triton* other) { double SE = r*r; return NULL; }
	virtual Products* products(Helion3* other) { double SE = 8.7e-41; return NULL; }
	virtual Products* products(Helion4* other) { return NULL; }

};

//Double-dispatch for Particle
void Particle::interact(Particle* other) { other->interact(this); }
void Particle::interact(Photon* other) { other->interact(this); }
void Particle::interact(Positron* other) { other->interact(this); }
void Particle::interact(Neutrino* other) { other->interact(this); }
void Particle::interact(Nucleon* other) { other->interact(this); }
void Photon::interact(Particle* other) { other->interact(this); }
void Positron::interact(Particle* other) { other->interact(this); }
void Neutrino::interact(Particle* other) { other->interact(this); }
void Nucleon::interact(Particle* other) { other->interact(this); }

//Double-dispatch for Nucleon::SE()
Products* Nucleon::products(Nucleon* other) { return other->products(this); }
Products* Nucleon::products(Proton* other) { return other->products(this); }
Products* Nucleon::products(Deuteron* other) { return other->products(this); }
Products* Nucleon::products(Triton* other) { return other->products(this); }
Products* Nucleon::products(Helion3* other) { return other->products(this); }

///Massive-massive interaction
//Elastic collision, bremsstrahlung radiation, and nucleosynthesis
void Nucleon::interact(Nucleon* other) {
	printf("%s interacted with %s\n", other->name,this->name);

	//Elastic collision
	double time_to_collide = abs(this->x - other->x) / get_v();
	int num_of_collisions = timestep / time_to_collide;

}

//Elastic collision and bremsstrahlung radiation
void Positron::interact(Positron* other) {
	printf("Positron interacted with Positron\n");
}

//Elastic collision and bremsstrahlung radiation
void Positron::interact(Nucleon* other) { other->interact(this); }
void Nucleon::interact(Positron* other) {
	printf("%s interacted with %s\n", other->name,this->name);
}

///Photonic interactions
//Two-photon pair-production
void Photon::interact(Photon* other) {

}

//Scattering (1D Klein-Nishina) and pair-production
void Positron::interact(Photon* other) { other->interact(this); }
void Photon::interact(Positron* other) {
	printf("%s interacted with %s\n", other->name,this->name);
}

//Scattering (1D Klein-Nishina) and pair-production
void Photon::interact(Nucleon* other) { other->interact(this); }
void Nucleon::interact(Photon* other) {
	printf("%s interacted with %s\n", other->name,this->name);
}

///Neutrino interactions
void Neutrino::interact(Neutrino* other) {} //Lorite creation?

void Positron::interact(Neutrino* other) { other->interact(this); }
void Neutrino::interact(Positron* other) {}

void Photon::interact(Neutrino* other) { other->interact(this); }
void Neutrino::interact(Photon* other) {}

void Nucleon::interact(Neutrino* other) { other->interact(this); }
void Neutrino::interact(Nucleon* other) {}

int main() {
	Particle* middle = new Photon(10 * e, 0);
	Particle* iter = middle;
	iter = iter->right = new Deuteron(10 * e, 1e-9);
	iter = iter->right = new Triton(10 * e, 2e-9);
    Particle* abcd = iter;
	iter = iter->right = new Neutrino(10 * e, 3e-9);
	iter = iter->right = new Proton(10 * e, 4e-9);

    Particles* particles = new Particles(middle);

    abcd->remove();
    printf("THIS %s\n", abcd->name);

    if (!particles->check_validity()) {
        printf("ERROR: Particles invalid");
        return 1;
    }

    particles->print();

    particles->interact();

	return 0;
}
