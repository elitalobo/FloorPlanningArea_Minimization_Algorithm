#include<bits/stdc++.h>
using namespace std;
#define XMAX 100
#define XMIN 0
#define	YMAX 100
#define YMIN 0
#define wt 0.01
#define Tv 0.1
#define inf 10000
#define ll long long
#define Rand() ((double)rand()/RAND_MAX)

struct Particle {
	double *x;
	double *v;
	double f;
	double pbest;
	double *x_star;
} *particles;

#define Nparticles 50
#define T_MAX 100000
#define W_0 0.9
#define W_T 0.4
#define MAX_V 5.0
#define c1 2.0
#define c2 2.0
// MAX_V was 2.0
#define Nvariables 2

#define better(y1,y2) (y1<y2)
int N;
struct Block {
        double x, y;
        double height, width;
	double pbest_x, pbest_y;
	double gbest_x, gbest_y;
        double area;
	double vel_x;
	double vel_y;
	Block() {}
	Block(int h, int w) {
		pbest_x = -inf;
		pbest_y = -inf;
		gbest_x = -inf;
		gbest_y = -inf;
		area = h*w;
		height = h;
		width=w;
		x=inf;
		y=inf;
		vel_x=inf;
		vel_y =inf;
	}
		
} blocks[1000];

bool compare(struct Block const &a, struct Block const &b) {
	if(a.area==b.area) {
		if(a.height==b.height)
			return a.width > b.width;
		return a.height > b.height;
	}
	return a.area > b.area;
}

void initialize_blocks(int n){
	for(int i=0;i<n;i++) {
		blocks[i].x = inf;
		blocks[i].y = inf;
		blocks[i].vel_x = 0.0;
		blocks[i].vel_y = 0.0;
		blocks[i].pbest_x = inf;
		blocks[i].pbest_x = inf;
	}
	
}



double find_area(int idx) {
	double max_x=0;
	double max_y =0;
	double min_x = inf;
	double min_y = inf;
	for(int i=0;i<idx;i++) {
		if(max_x < blocks[i].pbest_x) {
			max_x = blocks[idx].pbest_x;
		}
		if(min_x > blocks[i].pbest_x) {
			min_x = blocks[idx].pbest_x;
		}
		if(max_y < blocks[i].pbest_y) {
			max_y = blocks[idx].pbest_y;
		}
		if(min_y > blocks[i].pbest_y) {
			min_y = blocks[idx].pbest_y;
		}
	}
	return abs((max_x-min_x)*(max_y-min_y));
}

	
void print_best(int n) {
		for(int i=0;i<n;i++) {
			cout <<"block "  << i << "  should be placed at (" << blocks[i].pbest_x << "," << blocks[i].pbest_y << ") " << endl;
		}
		return;
}

void TakeInput() {
	scanf("%d",&N);
	for(int i=0;i<N;i++) {
		double x,y;
		scanf("%lf %lf",&x,&y);
		blocks[i] = Block(x,y);
	}
	initialize_blocks(N);
	sort(blocks,blocks+N,compare);
}

bool is_not_overlapping(Particle p, int idx) {
	double x_cord = p.x[0];
	double y_cord = p.x[1];
        for(int i=0;i<idx;i++) {
                bool flag=true;
                if(!((blocks[i].pbest_x + blocks[i].width <= x_cord) || (blocks[i].pbest_y + blocks[i].height <= y_cord) || (blocks[i].pbest_x - blocks[idx].width >= x_cord)|| (blocks[i].pbest_y - blocks[idx].height >= y_cord)))
                        {
                                return false;
                        }
        }
        return  true;
}

double find_area(Particle p,int idx) {
        double max_x=0;
        double max_y =0;
        double min_x = inf;
        double min_y = inf;
        for(int i=0;i<idx;i++) {
                if(max_x < blocks[i].pbest_x) {
                        max_x = blocks[idx].pbest_x;
                }
                if(min_x > blocks[i].pbest_x) {
                        min_x = blocks[idx].pbest_x;
                }
                if(max_y < blocks[i].pbest_y) {
                        max_y = blocks[idx].pbest_y;
                }
                if(min_y > blocks[i].pbest_y) {
                        min_y = blocks[idx].pbest_y;
                }
        }
	double x_cord = p.x[0];
	double y_cord = p.x[1];
	max_x = max(max_x,x_cord);
	max_y = max(max_y,y_cord);
	min_x = min(min_x, x_cord);
	min_y = min(min_y, y_cord);
        return abs((max_x-min_x)*(max_y-min_y));
}


void evaluate(Particle P, int idx) {
	P.f=0.0;
	if(!is_not_overlapping(P,idx))
	{
		return;
	}
	P.f = 1.0/find_area(P,idx);
}

void UpdateBest(Particle P) {
	for(int j=0;j<Nvariables;j++)
		P.x_star[j]=P.x[j];
	P.pbest = P.f;
}

int Initialize(int n, int l) {
	int G=0;
	for(int i=0;i<n;i++) {
		for(int j=0;j<Nvariables;j++) {
			particles[i].x[j]= XMIN + Rand()*(XMAX-XMIN);
			particles[i].v[j]= 0.0;
		}
		evaluate(particles[i],l);
		if(is_not_overlapping(particles[i],l)) {
			UpdateBest(particles[i]);
			if(better(particles[i].f,particles[G].f) && is_not_overlapping(particles[i],l))
				G=i;
		}
	}
	return G;
}


void NewParticles(int n) {
	int i;
	particles = (struct Particle*)malloc((n+1)*sizeof(struct Particle));
	for(int i=0;i<n;i++) {
		particles[i].x = (double*)malloc((Nvariables+1)*sizeof(double));
		particles[i].v = (double*)malloc((Nvariables+1)*sizeof(double));
		particles[i].x_star = (double*)malloc((Nvariables+1)*sizeof(double));
	}
	return;
}	
	

void Print(Particle P) {
	for(int i=0;i<Nvariables;i++)
		printf("%f",P.x_star[i]);
	printf(" =%lf\n",P.pbest);
}

int main() {
	int G;
	double w;
	w=W_0;
	TakeInput();
	blocks[0].pbest_x=0;
        blocks[0].pbest_y=0;
        blocks[0].x =0;
        blocks[0].y=0;
	for(int l=0;l<N;l++) {
		w=W_0;
		NewParticles(Nparticles);
		cout << "Done making new particles" << endl;
		Initialize(Nparticles,l);
		cout << "done initializing" << endl;
		G=0;
		for(int t=1;t<T_MAX;t++) {
			for(int i=0;i<Nparticles;i++) {
				for(int j=0;j<Nvariables;j++) {
					particles[i].v[j]=w*particles[i].v[j]+c1*Rand()*(particles[i].x_star[j]-particles[i].x[j])+c2*Rand()*(particles[G].x_star[j]-particles[i].x[j]);
				if(particles[i].v[j]<-MAX_V)
					particles[i].v[j]=-MAX_V;
				else if(particles[i].v[j]>MAX_V)
					particles[i].v[j]=MAX_V;
				particles[i].x[j]+=particles[i].v[j];
			}
			evaluate(particles[i],l);
			if(better(particles[i].f, particles[i].pbest) && is_not_overlapping(particles[i],l)) {
				if(better(particles[i].f, particles[G].pbest)) G=i;
					UpdateBest(particles[i]);
			}
		}
		//printf("%4d: ", t); Print(particles[G]);
		w-=(W_0-W_T)/T_MAX;
		}
		blocks[l].pbest_x = particles[G].x[0];
		blocks[l].pbest_y = particles[G].x[1];
	}
	print_best(N);
}



