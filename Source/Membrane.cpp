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
Membrane::Membrane (double T, double rho, double H, double kappaSq, double sig0, double sig1, double k, double Lx, double Ly) : T (T), rho (rho), H (H),
                  kappaSq (kappaSq),
                  sig0 (sig0),
                  sig1 (sig1),
                  Lx (Lx), Ly (Ly), k (k)
{
    excitationWidth = maxExcitationWidth;
    
    
    sig0Init = sig0;
    sig1Init = sig1;
    
    sig0 = 0;
    Tinit = T;
    cSq = T / (rho * H);
    h = 4.0 * 2.0 * sqrt ((cSq * k * k + 4.0 * sig1 * k + sqrt(pow(cSq * k * k + 4.0 * sig1 * k, 2) + 4.0 * kappaSq * k * k)) / 2.0);
    
//    h = sqrt((cSq * k * k + 4.0 * sig1 * k + sqrt(pow(cSq * k * k + 4.0 * sig1 * k, 2.0) + 16 * kappaSq * k * k)) / 2.0);
    
//    h = 0.062;
    Nx = floor (Lx/h);
    Ny = floor (Ly/h);
    h = Lx > Ly ? Lx / Nx : Ly / Ny;
    N = (Nx - 1) * (Ny - 1);
//    N = floor (aspectRatio / (h * h));
    
//    h = 1.0 / sqrt (N);
//    Nx = floor(sqrt (aspectRatio) / h);
//    Ny = floor(1.0 / (sqrt (aspectRatio) * h));
    
    
#ifdef ONEDVEC
//    for (int i = 0; i < numTimeSteps; ++i)
//        uVecs.push_back (std::vector<double> (N, 0));
//    uVecs.resize (numTimeSteps, std::vector<double> (N, 0)); //resize to three time steps
    uVecs.push_back (new std::vector<double> (N, 0));
    uNext = &uVecs[0][0][0];
    uVecs.push_back (new std::vector<double> (N, 0));
    u = &uVecs[1][0][0];
    uVecs.push_back (new std::vector<double> (N, 0));
    uPrev = &uVecs[2][0][0];
#else
    uVecs.resize (numTimeSteps, std::vector<std::vector<double>> (Nx, std::vector<double>(Ny, 0)));
    uNext = &uVecs[0][0];
    u = &uVecs[1][0];
    uPrev = &uVecs[2][0];
//    excited = true;
//    for (int i = 0; i < numTimeSteps; ++i)
//        uVecs.push_back (std::vector<std::vector<double>> (Nx, std::vector<double>(Ny, 0)));
#endif
    
    
    lambdaSq = cSq * k * k / (h * h);
    muSq = kappaSq * k * k / (h * h * h * h);
    
    d = 1.0 / (1.0 + sig0 * k);
    // u coeffs
    B1 = (2.0 - 4.0 * lambdaSq - 20.0 * muSq - 8.0 * sig1 * k / (h*h)) * d;
    B2 = (-muSq) * d;
    B3 = (-2.0 * muSq) * d;
    B4 = ((8.0 * muSq + lambdaSq + 2.0 * sig1 * k / (h*h))) * d;
    
    // uPrev coeffs
    C1 = (sig0 * k - 1 + 8.0 * sig1 * k / (h*h)) * d;
    C2 = (-2.0 * sig1 * k / (h*h)) * d;
    
    coeffs.push_back(B1);
    coeffs.push_back(B2);
    coeffs.push_back(B3);
    coeffs.push_back(B4);
    coeffs.push_back(C1);
    coeffs.push_back(C2);
    
//    d = 1.0f / (1.0f + sig0 * k);
//    A1 = cSq * k * k / (h * h) * d;
//    B1 = -(kappaSq * k * k) / (h * h * h * h) * d;
//    B2 = B1 * 2.0f;
//    B3 = B1 * -8.0f;
//    C = (2.0f * sig1 * k) / (h * h);
//    C1 = (2.0f - 4.0f * C + 20.0f * B1) * d;
//    C2 = (sig0 * k - 1.0f + 4.0f * C) * d;
//    C3 = C * d;
//    C4 = (k * k) * d;
//    excited = true;
    
}

Membrane::~Membrane()
{
    system("rm *.so");
    system("rm -rf *.so.dSYM");
    system("rm code.c");
}

void Membrane::paint (Graphics& g)
{
    float stateWidth = getWidth() / static_cast<double> (Nx - 4);
    float stateHeight = getHeight() / static_cast<double> (Ny - 4);
    int scaling = 10;
    for (int x = 2; x < Nx - 2; ++x)
    {
        for (int y = 2; y < Ny - 2; ++y)
        {
#ifdef ONEDVEC
            int cVal = clamp (255 * 0.5 * (u[x + y * Nx] * scaling + 1), 0, 255);
#else
            int cVal = clamp (255 * 0.5 * (u[x][y] * scaling + 1), 0, 255);
#endif
            g.setColour(Colour::fromRGB (cVal, cVal, cVal));
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
#ifdef ONEDVEC
    if (shouldBeZero)
    {
        for (int i = 0; i < N; ++i)
        {
            if (uNext[i] != 0)
            {
                std::cout << "uNext is not 0 at" + String(i) + "!" << std::endl;
            }
            if (u[i] != 0)
            {
                std::cout << "u is not 0 at" + String(i) + "!" << std::endl;
            }
            if (uPrev[i] != 0)
            {
                std::cout << "uPrev is not 0 at" + String(i) + "!" << std::endl;
            }
        }
    }
    updateEq (uNext, u, uPrev, &coeffs[0], Nx, Ny);
    if (shouldBeZero)
    {
        for (int i = 0; i < N; ++i)
        {
            if (uNext[i] != 0)
            {
                std::cout << "uNext is not 0 at" + String(i) + "!" << std::endl;
            }
            if (u[i] != 0)
            {
                std::cout << "u is not 0 at" + String(i) + "!" << std::endl;
            }
            if (uPrev[i] != 0)
            {
                std::cout << "uPrev is not 0 at" + String(i) + "!" << std::endl;
            }
        }
    }
//    if (!update)
//    {
//    for (int t = 0; t < uVecs.size(); ++t)
//    {
//        for (int i = 0; i < uVecs[0].size(); ++i)
//        {
//            if (uVecs[t][i] != 0)
//            {
//                std::cout << "Uvecs " + String(t) + " is not 0!" << std::endl;
//            }
//        }
//    }
//    else {
//        std::cout << "wait" << std::endl;
//    }
#else
    for (int l = 2; l < Nx - 2; ++l)
    {
        for (int m = 2; m < Ny - 2; ++m)
        {
            uNext[l][m] = B1 * u[l][m]
            + B2 * (u[l+2][m] + u[l-2][m] + u[l][m+2] +  u[l][m-2])
            + B3 * (u[l+1][m+1] + u[l+1][m-1] + u[l-1][m+1] +  u[l-1][m-1])
            + B4 * (u[l+1][m] + u[l-1][m] + u[l][m+1] +  u[l][m-1])
            + C1 * uPrev[l][m]
            + C2 * (uPrev[l+1][m] + uPrev[l-1][m] + uPrev[l][m+1] + uPrev[l][m-1]);
        }
    }
#endif
}

void Membrane::updateStates()
{
//    int lengthUVec = static_cast<int>(u.size());
//    for (int i = lengthUVec - 1; i > 0; --i)
//        u[i] = u[i - 1];
//    uNextPtrIdx = (uNextPtrIdx + (lengthUVec - 1)) % lengthUVec;
//    u[0] = &uVecs[uNextPtrIdx][0];
    
    if (update)
    {
#ifdef ONEDVEC
        double* dummyPtr = uPrev;
#else
        std::vector<double>* dummyPtr = uPrev;
#endif
        uPrev = u;
        u = uNext;
        uNext = dummyPtr;
    }
}


void Membrane::excite()
{
    excitationWidth = clamp (excitationWidth, 1, Nx - 5);
    if (excited)
    {
        shouldBeZero = false;
        excited = false;
        excitationTime = Time::currentTimeMillis();

        std::vector<std::vector<double>> excitationArea (excitationWidth, std::vector<double> (excitationWidth, 0));
        
        for (int i = 1; i < excitationWidth; ++i)
        {
            for (int j = 1; j < excitationWidth; ++j)
            {
                excitationArea[i][j] = 10.0 / (excitationWidth * excitationWidth) * 0.25 * (1 - cos(2.0 * double_Pi * i / static_cast<int>(excitationWidth+1))) * (1 - cos(2.0 * double_Pi * j / static_cast<int>(excitationWidth+1)));
            }
        }
        
        int startIdX = clamp(idX - excitationWidth * 0.5, 2, Nx-3 - excitationWidth);
        int startIdY = clamp(idY - excitationWidth * 0.5, 2, Ny-3 - excitationWidth);
//
//        for (int i = 1; i < excitationWidth; ++i)
//        {
//            for (int j = 1; j < excitationWidth; ++j)
//            {
//#ifdef ONEDVEC
//                u[static_cast<int> (i + startIdX + Nx * (j + startIdY))] += excitationArea[i][j];
//                uPrev[static_cast<int> (i + startIdX + Nx * (j + startIdY))] += excitationArea[i][j];
//#else
//                u[i + startIdX][j + startIdY] += excitationArea[i][j];
//                uPrev[i + startIdX][j + startIdY] += excitationArea[i][j];
                u[7][7] += 1.0;
                uPrev[7][7] += 1.0;

//#endif
//            }
//        }
//        uPrev[static_cast<int> (clamp(idX, 2, Nx-3))][static_cast<int> (clamp(idY, 2, Ny-3))] += 1.0;
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
    idX = e.x / stateWidth + 2;
    idY = e.y / stateHeight + 2;
    idX = 7;
    idY = 7;
    
}
void Membrane::mouseDrag (const MouseEvent& e)
{
    int stateWidth = getWidth() / static_cast<double> (Nx - 4);
    int stateHeight = getHeight() / static_cast<double> (Ny - 4);
    idX = e.x / static_cast<double>(stateWidth) + 2;
    idY = e.y / static_cast<double>(stateHeight) + 2;
    idX = 7;
    idY = 7;
    excited = true;
}

void Membrane::mouseUp (const MouseEvent &e)
{
}

void Membrane::createUpdateEq()
{
    void *handle;
    char *error;

    // convert updateEqString to char
    
    std::hash<int64> hasher;
    auto newName = hasher (juce::Time::getCurrentTime().toMilliseconds());
    
    FILE *fd= fopen("code.c", "w");
    
    const char* eq =
    " for (int l = 2; l < Nx - 2; ++l)\n"
    "{ \n"
    "for (int m = 2; m < Ny - 2; ++m)\n"
    "{ \n"
    "uNext[l + m * Nx] ="
    "   coeffs[0] * u[l + m * Nx]"
    " + coeffs[1] * (u[l+2 + m * Nx] + u[l-2 + m * Nx] + u[l + (m+2) * Nx] + u[l + (m-2) * Nx])"
    " + coeffs[2] * (u[l+1 + (m+1) * Nx] + u[l-1 + (m-1) * Nx] + u[l+1 + (m-1) * Nx] + u[l-1 + (m+1) * Nx])"
    " + coeffs[3] * (u[l+1 + m * Nx] + u[l-1 + m * Nx] + u[l + (m+1) * Nx] + u[l + (m-1) * Nx])"
    " + coeffs[4] * uPrev[l + m * Nx]"
    " + coeffs[5] * (uPrev[l+1 + m * Nx] + uPrev[l-1 + m * Nx] + uPrev[l + (m+1) * Nx] + uPrev[l + (m-1) * Nx]);\n"
    "} \n"
    "}";
    
    fprintf(fd, "#include <stdio.h>\n"
            "void updateEq(double* uNext, double* u, double* uPrev, double* coeffs, int Nx, int Ny)\n"
            "{\n"
            "%s\n"
            "}", eq);
    fclose(fd);
    
    String systemInstr;
    
    // if this object is new (has not been edited yet) delete the previous file
    systemInstr = String ("clang -shared -undefined dynamic_lookup -O3 -o " + String (newName) + ".so code.c -g");
    system (toConstChar (systemInstr));
    handle = dlopen (toConstChar (String (String (newName) + ".so")), RTLD_LAZY);
   
    if (!handle)
    {
        fprintf (stderr, "%s\n", dlerror());
        exit(1);
    }
    
    dlerror();    /* Clear any existing error */
    
    *(void **)(&updateEq) = dlsym (handle, "updateEq"); // second argument finds function name
    
    if ((error = dlerror()) != NULL)  {
        fprintf (stderr, "%s\n", error);
        exit(1);
    }
}

void Membrane::setTension (double newTension)
{
    T = newTension;
    cSq = T / (rho * H);
    lambdaSq = cSq * k * k / (h * h);
    updateParams();
#ifdef ONEDVEC
    coeffs[0] = B1;
    coeffs[3] = B4;
#endif
}

void Membrane::setSig0 (double newSig0)
{
    sig0 = newSig0;
    updateParams();
}

void Membrane::setSig1 (double newSig1)
{
    sig1 = newSig1;
    updateParams();
}

void Membrane::updateParams()
{
    d = 1.0 / (1.0 + sig0*k);
    // u coeffs
    B1 = (2.0 - 4.0 * lambdaSq - 20.0 * muSq - 8.0 * sig1 * k / (h*h)) * d;
    B2 = (-muSq) * d;
    B3 = (-2.0 * muSq) * d;
    B4 = ((8.0 * muSq + lambdaSq + 2.0 * sig1 * k / (h*h))) * d;
    
    // uPrev coeffs
    C1 = (sig0 * k - 1 + 8.0 * sig1 * k / (h*h)) * d;
    C2 = (-2.0 * sig1 * k / (h*h)) * d;
    
    coeffs[0] = B1;
    coeffs[1] = B2;
    coeffs[2] = B3;
    coeffs[3] = B4;
    coeffs[4] = C1;
    coeffs[5] = C2;
}

void Membrane::changeDamping()
{
    sig1 = sig1Init * exp(-alpha * (Time::currentTimeMillis() - excitationTime));
    updateParams();
    
}
