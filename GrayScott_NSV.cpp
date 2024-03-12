#include "Declarations.h"
#include "Headers/NSVClasses.h"
#include "Headers/GlobalParameters.h"

using namespace std;

// to generate random numbers
namespace myrand
{
  std::mt19937 rng;
}

void GrayScottSimulation(vector<ofstream> *outputs)
{
    // initialization
    // I have 2 morphogens (species): U and V
    vector<double>u(N,ceil(1.0*Omega*h));
    vector<double>v(N,ceil(1.0*Omega*h));
    Propensities propensities(R,N);
    EventQueue events;
    Transition currentTransition;

    int voxel=0, destination_voxel=0, reactionID=0;
    bool jump_event=0;
    int count_print=0, i=0, count=0;
    double time=0.0, timeprint = Tmax/(double)200.0, tau=0.0;

    // calculate initial propensities (this creates initial event queue)
    for(i=0;i<N;i++)
        UpdatePropensities(&u,&v,&propensities,&events,i+1,false,time);

    // beginning of the main loop (just as in Gillespie)
    while(time<Tmax)
    {
        count++;
        //print out into the files
        if(time>count_print*timeprint)
        {
            count_print++;
            print_vector(&u, &outputs->at(0), time, true);
            print_vector(&v, &outputs->at(1), time, true);
            cout<<time<<endl;
        }

        // get next transition (in which voxel and at what time)
        currentTransition=events.top();
        events.pop();
        voxel=currentTransition.id;
        tau=currentTransition.time-time;
        time+=tau;
        if(time>Tmax)
            goto EndOfSimulation;
        // this probabilistically decides which reactions occurs in the given voxel
        reactionID=propensities.wchEvent(voxel);
        
        // here you should implement the corresponding changes in your species for each reaction
        switch(reactionID)
        {
            case 1:     // U + 2V ->[rho] 3V
            {
                u[voxel-1]--;
                v[voxel-1]++;
                jump_event=0;
                break;
            }
            case 2:     // 0 ->[F] U
            {
                u[voxel-1]++;
                jump_event=0;
                break;
            }
            case 3:     // U ->[F] 0
            {
                u[voxel-1]--;
                jump_event=0;
                break;
            }
            case 4:     // V ->[F+k] 0
            {
                v[voxel-1]--;
                jump_event=0;
                break;
            }
            case 5:     // diffusion of U
            {
                u[voxel-1]--;
                jump_event=1;
                if(DoubleUniformDistribution(0.0,1.0)<0.5)  // jump to right
                {
                    if(voxel<N)
                    {
                        u[voxel]++;
                        destination_voxel=voxel+1;
                    }
                    else    // voxel=N
                    {
                        u[0]++;     // (periodic boundary condition)
                        destination_voxel=1;
                    }
                }
                else // jump to left
                {
                    if(voxel>1)
                    {
                        u[voxel-2]++;
                        destination_voxel=voxel-1;
                    }
                    else    // voxel=1
                    {
                        u[N-1]++;     // (periodic boundary condition)
                        destination_voxel=N;
                    }
                }
                break;
            }
            case 6:     // diffusion of V
            {
                v[voxel-1]--;
                jump_event=1;
                if(DoubleUniformDistribution(0.0,1.0)<0.5)  // jump to right
                {
                    if(voxel<N)
                    {
                        v[voxel]++;
                        destination_voxel=voxel+1;
                    }
                    else    // voxel=N
                    {
                        v[0]++;     // (periodic boundary condition)
                        destination_voxel=1;
                    }
                }
                else // jump to left
                {
                    if(voxel>1)
                    {
                        v[voxel-2]++;
                        destination_voxel=voxel-1;
                    }
                    else    // voxel=1
                    {
                        v[N-1]++;     // (periodic boundary condition)
                        destination_voxel=N;
                    }
                }
                break;
            }
            default:
            {
                cout<<"Error in reactionID"<<endl;
                goto EndOfSimulation;
                break;
            }
        }

        // update after each iteration only in the voxels which changed
        UpdatePropensities(&u,&v,&propensities,&events,voxel,false,time);
        if(jump_event)
            UpdatePropensities(&u,&v,&propensities,&events,destination_voxel,true,time);
    }

EndOfSimulation:
    // final printing out into the files
    print_vector(&u, &outputs->at(0), time, true);
    print_vector(&v, &outputs->at(1), time, true);
    for(int o=0;o<outputs->size();o++)
		(*outputs)[o].flush();
    return;
}

void UpdatePropensities(vector<double> *u, vector<double> *v,Propensities *propensities, EventQueue *events, int voxel, bool erase, double time)
{
    // if the voxel exists already, you have to delete it first from the event queue
	if(erase)
		events->erase(voxel);

    // Here you put your reactions:
    // you can make it more efficient by recalculating only the terms which have changed

    propensities->Reactions[0][voxel-1] = rho/Omega/Omega/h/h*(*u)[voxel-1]*(*v)[voxel-1]*((*v)[voxel-1]-1);  // U + 2V ->[rho] 3V
    propensities->Reactions[1][voxel-1] = F*Omega*h;            // 0 ->[F] U 
    propensities->Reactions[2][voxel-1] = F*(*u)[voxel-1];      // U ->[F] 0
    propensities->Reactions[3][voxel-1] = (F+k)*(*v)[voxel-1];  // V ->[F+k] 0
    propensities->Reactions[4][voxel-1] = 2.0*Du/h/h*(*u)[voxel-1]; // diffusion of U
    propensities->Reactions[5][voxel-1] = 2.0*Dv/h/h*(*v)[voxel-1]; // diffusion of V

    // you do not need to modify the code below
    // it calculates the total propensity in the voxel and inserts it into the event queue
    double TotalPropensity = 0.0, tau=0.0, rand=DoubleUniformDistribution(0, 1);

    for(int i=0; i<propensities->R; i++)
		TotalPropensity+=propensities->Reactions[i][voxel-1];
	
	if (TotalPropensity > 0.0)
	{
		tau = ((double)1 / (double)TotalPropensity )*log(1 / rand);
		events->push(Transition(voxel, time + tau));
	}
}

// for printing out vectors into a file
void print_vector(vector<double> *v, ofstream *output, double time, double print_time)
{
	copy(v->begin(), v->end(), ostream_iterator<double>(*output, " "));
    if(print_time)
        *output<<time<<endl;
    else
        *output<<endl;
}

// to generate a uniformly distributed random number in the interval [min,max]
double DoubleUniformDistribution(double min, double max)
{
	return min + (max - min) * (((double)rand()) / ((double)RAND_MAX));
}
