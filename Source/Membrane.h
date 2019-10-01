#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>
#include <string.h>
//==============================================================================
/*
 */
#define ONEDVEC

class Membrane    : public Component
{
public:
    Membrane (double cSq, double kappaSq, double sig0, double sig1, double k, double Lx, double Ly);
    ~Membrane();
    
    void paint (Graphics& g) override;
    void resized() override;
    
    void calculateFDS();
#ifdef ONEDVEC
    double getOutput (double xRatio, double yRatio) {return u[static_cast<int> (xRatio * Nx + Nx * yRatio * Ny)]; };
#else
    double getOutput (double xRatio, double yRatio) {return u[static_cast<int> (xRatio * Nx)][static_cast<int>(yRatio * Ny)]; };
#endif
    void updateStates();
    
    void setImpactPosition (float xPos, float yPos);
    void excite();
//    void setZero();
    
    void setAspectRatio (double aspRat) { aspectRatio = aspRat; }; //probably to some state-interpolation or refreshing...
    
    double clamp (double val, double min, double max);
    
    void mouseDown (const MouseEvent& e) override;
    void mouseDrag (const MouseEvent& e) override;
    void mouseUp (const MouseEvent& e) override;
    
    bool isExcited() { return excited; };
    
    void createUpdateEq();
private:
    double cSq, kappaSq, rho, E, nu, H, sig0, sig1, Lx, Ly; // material values
    double lambdaSq, muSq; // courant stuffs
    double d, A1, B1, B2, B3, B4, C, C1, C2;
    
    double k; // timestep
    
    // pointers to the different states
#ifdef ONEDVEC
    double* uNext;
    double* u;
    double* uPrev;
#else
    std::vector<double>* uNext;
    std::vector<double>* u;
    std::vector<double>* uPrev;
#endif
    std::vector<double> coeffs;
    // states
#ifdef ONEDVEC
    std::vector<std::vector<double>> uVecs;
#else
    std::vector<std::vector<std::vector<double>>> uVecs;
#endif
    
    
    double aspectRatio = 1;
    int Nx;
    int Ny;
    
    int N;
    double h;
    
    int numTimeSteps = 3;
    
    bool excited = false;
    
    int idX, idY;
    void (*updateEq) (double* uNext, double* u, double* uPrev, double* parameters, int Nx, int Ny);
    

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (Membrane)
};
