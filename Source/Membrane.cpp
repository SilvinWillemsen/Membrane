/*
 ==============================================================================
 
 Membrane.cpp
 Created: 24 Jul 2019 3:48:59pm
 Author:  Silvin Willemsen
 
 ==============================================================================
 */

#include "../JuceLibraryCode/JuceHeader.h"
#include "Membrane.h"

//==============================================================================
Membrane::Membrane (double cSq, double kappaSq, double sig0, double sig1, double k) : cSq (cSq),
                                                                                      kappaSq (kappaSq),
                                                                                      sig0 (sig0),
                                                                                      sig1 (sig1),
                                                                                      k (k)
{
    
    h = sqrt((cSq * k * k + 4.0 * sig1 * k + sqrt(pow(cSq * k * k + 4.0 * sig1 * k, 2.0) + 16 * kappaSq * k * k)) / 2.0);
    
    h = 0.062;
    N = floor (aspectRatio / (h * h));
    
    h = 1.0 / sqrt (N);
    Nx = floor(sqrt (aspectRatio) / h);
    Ny = floor(1.0 / (sqrt (aspectRatio) * h));
    
    uVecs.reserve (numTimeSteps); //resize to three time steps
    
    for (int i = 0; i < numTimeSteps; ++i)
        uVecs.push_back (std::vector<std::vector<double>> (Ny, std::vector<double>(Nx, 0)));
    
    
    uNext = &uVecs[0][0];
    u = &uVecs[1][0];
    uPrev = &uVecs[2][0];
    
    
    muSq = kappaSq * k * k / (h * h * h * h);
    
    d = 1.0f / (1.0f + sig0 * k);
    A1 = cSq * k * k / (h * h) * d;
    B1 = -(kappaSq * k * k) / (h * h * h * h) * d;
    B2 = B1 * 2.0f;
    B3 = B1 * -8.0f;
    C = (2.0f * sig1 * k) / (h * h);
    C1 = (2.0f - 4.0f * C + 20.0f * B1) * d;
    C2 = (sig0 * k - 1.0f + 4.0f * C) * d;
    C3 = C * d;
    C4 = (k * k) * d;
    
}

Membrane::~Membrane()
{
}

void Membrane::paint (Graphics& g)
{
    float stateWidth = getWidth() / static_cast<double> (Nx - 4);
    float stateHeight = getHeight() / static_cast<double> (Ny - 4);
    int scaling = 10000;
    for (int x = 2; x < Nx - 2; ++x)
    {
        for (int y = 2; y < Ny - 2; ++y)
        {
            int cVal = clamp (255 * 0.5 * (u[x][y] * scaling + 1), 0, 255);
            g.setColour(Colour::fromRGBA (cVal, cVal, cVal, 127));
            g.fillRect((x - 2) * stateWidth, (y - 2) * stateHeight, stateWidth, stateHeight);
        }
    }
}

void Membrane::resized()
{
}

void Membrane::calculateFDS()
{
    // input states and coefficients in vector form
    for (int m = 2; m < Ny - 2; ++m)
    {
        for (int l = 2; l < Nx - 2; ++l)
        {
            uNext[l][m] = A1 * (u[l][m+1] + u[l][m-1] + u[l+1][m] + u[l-1][m] - 4.0 * u[l][m]) +
            B1 * (u[l][m+2] + u[l][m-2] + u[l+2][m] + u[l-2][m])
            + B2 * (u[l+1][m+1] + u[l-1][m+1] + u[l+1][m-1] + u[l-1][m-1])
            + B3 * (u[l][m+1] + u[l][m-1] + u[l+1][m] + u[l-1][m])
            + C1 * u[l][m]
            + C2 * uPrev[l][m];
            //                    + C3 * (u[l + (m+1) * Nx] + u[l + (m-1) * Nx] + u[l+1 + m * Nx] + u[l-1 + m * Nx] - uPrev[l + (m+1) * Nx] - uPrev[l + (m-1) * Nx] - uPrev[l+1 + m * Nx] - uPrev[l-1 + m * Nx]);
            //            std::cout << l + m * Nx << std::endl;
        }
    }
    
}

void Membrane::updateStates()
{
//    int lengthUVec = static_cast<int>(u.size());
//    for (int i = lengthUVec - 1; i > 0; --i)
//        u[i] = u[i - 1];
//    uNextPtrIdx = (uNextPtrIdx + (lengthUVec - 1)) % lengthUVec;
//    u[0] = &uVecs[uNextPtrIdx][0];
    
    std::vector<double>* dummyPtr = uPrev;
    uPrev = u;
    u = uNext;
    uNext = dummyPtr;
}


void Membrane::excite()
{
    if (excited)
    {
        excited = false;
        u[static_cast<int>(0.5*Nx)][static_cast<int>(0.5 * Ny)] = 1.0;
        uPrev[static_cast<int>(0.5*Nx)][static_cast<int>(0.5*Ny)] = 1.0;
        
//        int width = floor ((N * 2.0) / 5.0) / 2.0;
//        int loc = floor (N * static_cast<float>(getXLoc()) / static_cast<float>(getWidth()));
//        int startIdx = clamp (loc - width / 2.0, simplySupported ? 1 : 2, simplySupported ? N-1-width : N-2-width);
//        for (int i = 0; i < width; ++i)
//        {
//            double val = (1 - cos (2 * double_Pi * i / width)) * 0.5;
//            for (int j = 1; j < u.size(); ++j)
//                u[j][startIdx + i] = u[j][startIdx + i] + val;
//        }
    }
}

double Membrane::clamp (double val, double min, double max)
{
    if (val < min)
    {
        val = min;
        return val;
    }
    else if (val > max)
    {
        val = max;
        return val;
    }
    return val;
}

//void Membrane::setImpactPosition (float xPos, float yPos)
//{
//    double pointX = xPos * Nx;
//    double pointY = yPos * Ny;
//    strikePositionX = clamp(pointX, 2, Nx - 2);
//    strikePositionY = clamp(pointY, 2, Ny - 2);
//
//    int spX = floor (strikePositionX);
//    int spY = floor (strikePositionY);
//
//    for (int x = 0; x < Nx; x++)
//    {
//        for (int y = 0; y < Ny; y++)
//        {
//            if (interpolation == noPlateInterpol)
//            {
//                if (x == spX && y == spY)
//                    excitationArea[x][y] = 1.0f;
//                else
//                    excitationArea[x][y] = 0.0f;
//            } else if (interpolation == bilinear)
//            {
//                double alphaX = strikePositionX - spX;
//                double alphaY = strikePositionY - spY;
//
//                if (x == spX && y == spY)
//                {
//                    excitationArea[x][y] = (1-alphaX) * (1-alphaY) * 1.0f;
//                    if (spX < Nx - 3)
//                        excitationArea[x+1][y] = alphaX * (1-alphaY) * 1.0f;
//                    if (spY < Ny - 3)
//                        excitationArea[x][y+1] = (1-alphaX) * alphaY * 1.0f;
//                    if (spX < Nx - 3 && spY < Ny - 3)
//                        excitationArea[x+1][y+1] = alphaX * alphaY * 1.0f;
//                }
//            }
//        }
//    }
//}

void Membrane::mouseDown (const MouseEvent &e)
{
    int stateWidth = getWidth() / static_cast<double> (Nx - 4);
    int stateHeight = getHeight() / static_cast<double> (Ny - 4);
    int idX = e.x / stateWidth + 2;
    int idY = e.y / stateHeight + 2;
    
}
void Membrane::mouseDrag (const MouseEvent& e)
{
    int stateWidth = getWidth() / static_cast<double> (Nx - 4);
    int stateHeight = getHeight() / static_cast<double> (Ny - 4);
    double idX = e.x / static_cast<double>(stateWidth) + 2;
    double idY = e.y / static_cast<double>(stateHeight) + 2;
    
    excited = true;
}

void Membrane::mouseUp (const MouseEvent &e)
{
}

