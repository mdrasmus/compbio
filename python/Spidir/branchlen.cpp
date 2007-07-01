/*=============================================================================

    SPIDIR    
    Branch Length Estimation
    
    branchlen.cpp
    started: Sun Jul  1 13:11:02 EDT 2007

=============================================================================*/



#include "common.h"



/*=============================================================================
    From: Felsenstein. Inferring Phylogenies. p 202.

    NOTE: HKY is Tamura-Nei where 
        alpha_r / alpha_y = rho = pi_r / pi_y 

    NOTE: parameters are chosen such that 
        P(i != j | j, t=1, pi, ratio) = 1
        thus, time also corresponds to sub/site

    Definitions:
        i    = destination base
        j    = source base
        t    = time
        pi_i = background/prior distribution of base i
        R    = Transtition/Transversion ratio

    Parameterization:
        beta = 1 / (2.0 * pi_r * pi_y * (1+R))
        rho = pi_r / pi_y    

        alpha_y = (pi_r * pi_y * R - pi_a*pi_g - pi_c*pi_t) / 
                  (2.0*(1+R)*(pi_y*pi_a*pi_g*rho + pi_r*pi_c*pi_t))
        alpha_r = rho * alpha_y    

    Convenience variables:
        if (dnatype[i] == DNA_PURINE) {
            alpha_i = alpha_r;
            pi_ry = pi_r;
        } else {
            alpha_i = alpha_y;
            pi_ry = pi_y;
        }
        int delta_ij =  int(i == j);
        int e_ij = int(dnatype[i] == dnatype[j]);    

    prob(i | j, t, pi, R) = 
        exp(-(alpha_i + beta)*t) * delta_ij + 
        exp(-beta*t) * (1 - exp(-alpha_i*t)) * (pi_j*e_ij/pi_ry) + 
        (1 - exp(-beta*t)) * pi_j

*/
class HkyModel
{
public:
    HkyModel(float *bgfreq, float ratio) :
        ratio(ratio)
    {
        // set background base frequencies
        for (int i=0; i<4; i++)
            pi[i] = bgfreq[i];
        
        pi_r = pi[DNA_A] + pi[DNA_G];
        pi_y = pi[DNA_C] + pi[DNA_T];
        rho = pi_r / pi_y;


        // determine HKY parameters alpha_r, alpha_y, and beta
        b = 1.0 / (2.0 * pi_r * pi_y * (1.0+ratio));
        a_y = (pi_r*pi_y*ratio - pi[DNA_A]*pi[DNA_G] - pi[DNA_C]*pi[DNA_T]) / 
          (2.0*(1+ratio)*(pi_y*pi[DNA_A]*pi[DNA_G]*rho + pi_r*pi[DNA_C]*pi[DNA_T]));
        a_r = rho * a_y;
        
    }
    

    // transition probability
    float operator()(int i, int j, float t)
    {
        float pi_j = pi[j];
        
        //  convenience variables
        float a_i, pi_ry;
        if (dnatype[i] == DNA_PURINE) {
            a_i = a_r;
            pi_ry = pi_r;
        } else {
            a_i = a_y;
            pi_ry = pi_y;
        }
        int delta_ij =  int(i == j);
        int e_ij = int(dnatype[i] == dnatype[j]);

        // transition probability
        float prob = exp(-(a_i + b)*t) * delta_ij + 
                     exp(-b*t) * (1 - exp(-a_i*t)) * (pi_j*e_ij/pi_ry) + 
                     (1 - exp(-b*t)) * pi_j;
        return prob;        
    }
    
    // parameters
    float ratio;
    float pi[4];
    float pi_r;
    float pi_y;
    float rho;
    float b;
    float a_y;
    float a_r;
};


template <class Model>
float distLikelihood(float **probs1, float **probs2, 
                      int seqlen, Model &model, float *bgfreq, float t)
{
    // trivial case
    if (t < 0)
        return 0;
    
    float logl = 0;
    
    // interate over sequence
    for (int k=0; k<seqlen; k++) {
        float prob = 0.0;
        
        // integrate over all transitions
        for (int i=0; i<4; i++) {
            float term = bgfreq[i] * probs1[k][i];
            for (int j=0; j<4; j++)
                prob += term * probs2[k][j] * model(i, j, t);
        }
        logl += log(prob);
    }
    
    return logl;
}



template <class Func>
float secantRoot(Func f, float x0, float x1, int maxiter, 
                 float minx=.0001, float esp=.0001)
{
    float f0 = f(x0);
    for (int i=0; i<maxiter; i++) {
        if (fabs(x1 - x0) < esp)
            return x0;
        float f1 = f(x1);
        float x2 = x1 - (x1 - x0) * f1 / (f1 - f0);
        
        x0 = (x1 < minx) ? minx : x1;
        x1 = (x2 < minx) ? minx : x2;
        f0 = f1;
    }

    return x1;
}


template <class Model>
class DistLikelihoodFunc
{
public:
    DistLikelihoodFunc(float **probs1, float **probs2, int seqlen, 
                       float *bgfreq, Model *model, float step=.01) :
        probs1(probs1),
        probs2(probs2),
        bgfreq(bgfreq),
        model(model),
        step(step)
    {}


    float operator()(float t)
    {
        float y0 = distLikelihood(probs1, probs2, seqlen, *model, bgfreq, t);
        float y1 = distLikelihood(probs1, probs2, seqlen, *model, bgfreq, t + step);
        return (y1 - y0) / step;
    }
    
    float **probs1;
    float **probs2;
    int seqlen;
    float *bgfreq;
    Model *model;
    float step;
};


template <class Model>
float mleDistance(float **probs1, float **probs2, int seqlen, 
                  float *bgfreq, Model model, 
                  float t0=.1, float t1=1, float step=.01)
{
    const int maxiter = 20;
    DistLikelihoodFunc<Model> df(probs1, probs2, seqlen, bgfreq, &model, step);
    return secantRoot(df, t0, t1, maxiter);
}





float f(float **probs1, float **probs2, int selen, float *bgfreq, float ratio,
      int seqlen, float t)
{
    HkyModel hky(bgfreq, ratio);
    
    return mleDistance(probs1, probs2, seqlen, bgfreq, hky);
}
