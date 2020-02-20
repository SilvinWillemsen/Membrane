#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
//==============================================================================
/*
 */
//#define ONEDVEC

class Membrane    : public Component
{
public:
    Membrane (double T, double rho, double H, double kappaSq, double sig0, double sig1, double k, double Lx, double Ly);
    ~Membrane();
    
    void paint (Graphics& g) override;
    void resized() override;
    
    void calculateFDS();
#ifdef ONEDVEC
    double getOutput (double xRatio, double yRatio) {return soundScalar * u[static_cast<int> (xRatio * Nx + Nx * yRatio * Ny)]; };
#else
    double getOutput (double xRatio, double yRatio) {return soundScalar * u[static_cast<int> (xRatio * Nx)][static_cast<int>(yRatio * Ny)]; };
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
    
    void setExcited(float amp = 1.0) { excited = true; excitationGain = amp; };
    bool isExcited() { return excited; };
        
    double getMaxT() { return Tinit; };
    double getMaxSize() { return LxInit; };
    double getMaxSig0() { return sig0Init; };
    double getMaxSig1() { return sig1Init; };
    int getMaxExcitationWidth () { return maxExcitationWidth; };
    
    void setTension (double newTension);
    void setTensionFromSize (double Lx);

    void setSig0 (double newSig0);
    void setSig1 (double newSig1);
    void setExcitationWidth (int newWidth) { excitationWidth = newWidth; };
    
    double getTension() { return T; };
    double getSig0() { return sig0; };
    double getSig1() { return sig1; };
    
    void changeDamping();
    void updateParams();
    
    const char* toConstChar (String string) { return static_cast<const char*> (string.toUTF8()); }
    
    void setUpdate (bool b) {
//        std::vector<std::vector<double>> emptyVec (3, std::vector<double>(N, 0));
//        uVecs = emptyVec;
        update = true;
        
    };
    void setState (double val) {
//        for (int t = 0; t < 3; ++t)
//        {
            for (int i = 2; i < Nx - 2; ++i)
            {
                for (int j = 2; j < Ny - 2; ++j)
                {
//                    uVecs[t][i][j] = val;
                    uNext[i][j] = val;
                    u[i][j] = val;
                    uPrev[i][j] = val;
                }
            }
//        }
        
    };
    
    void setExcitationGain (double val) { excitationGain = val; };
    
private:
    double Tinit, sig0Init, sig1Init, T, rho, H, LxInit;
    double cSq, kappaSq, E, nu, sig0, sig1, Lx, Ly; // material values
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
//    std::vector<std::vector<double>> uVecs {3, std::vector<double> (576, 0)};
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
    
    int idX = 0;
    int idY = 0;
    void (*updateEq) (double* uNext, double* u, double* uPrev, double* parameters, int Nx, int Ny);
    
    double soundScalar = 5.0;
    
    bool dynamicDamping = true;
    
    int64 excitationTime = 0;
    double alpha = 0.01;
    
    int maxExcitationWidth = 10;
    int excitationWidth;
    
    bool update = true;
    bool shouldBeZero = true;
    
    double excitationGain = 1.0;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (Membrane)
};
