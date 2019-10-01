/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
MembraneAudioProcessor::MembraneAudioProcessor()
#ifndef JucePlugin_PreferredChannelConfigurations
     : AudioProcessor (BusesProperties()
                     #if ! JucePlugin_IsMidiEffect
                      #if ! JucePlugin_IsSynth
                       .withInput  ("Input",  AudioChannelSet::stereo(), true)
                      #endif
                       .withOutput ("Output", AudioChannelSet::stereo(), true)
                     #endif
                       )
#endif
{
}

MembraneAudioProcessor::~MembraneAudioProcessor()
{
}

//==============================================================================
const String MembraneAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool MembraneAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool MembraneAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

bool MembraneAudioProcessor::isMidiEffect() const
{
   #if JucePlugin_IsMidiEffect
    return true;
   #else
    return false;
   #endif
}

double MembraneAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int MembraneAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int MembraneAudioProcessor::getCurrentProgram()
{
    return 0;
}

void MembraneAudioProcessor::setCurrentProgram (int index)
{
}

const String MembraneAudioProcessor::getProgramName (int index)
{
    return {};
}

void MembraneAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void MembraneAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
    fs = sampleRate;
    double k = 1.0 / fs;
    double rho = 10;
    double E = 2e3;
    double H = 0.001;
    double nu = 0.3;
    double T = 100.0;
    double cSq = T / (rho * H);
    double kappaSq = E * H * H / (12.0 * rho * (1.0 - nu * nu));
    
    double sig0 = 0.1;
    double sig1 = 0.01;
    
    double Lx = 0.6;
    double Ly = 0.3;
    
    membranes.add (new Membrane (cSq, kappaSq, sig0, sig1, k, Lx, Ly));
    
}

void MembraneAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool MembraneAudioProcessor::isBusesLayoutSupported (const BusesLayout& layouts) const
{
  #if JucePlugin_IsMidiEffect
    ignoreUnused (layouts);
    return true;
  #else
    // This is the place where you check if the layout is supported.
    // In this template code we only support mono or stereo.
    if (layouts.getMainOutputChannelSet() != AudioChannelSet::mono()
     && layouts.getMainOutputChannelSet() != AudioChannelSet::stereo())
        return false;

    // This checks if the input layout matches the output layout
   #if ! JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;
   #endif

    return true;
  #endif
}
#endif

void MembraneAudioProcessor::processBlock (AudioBuffer<float>& buffer, MidiBuffer& midiMessages)
{
    ScopedNoDenormals noDenormals;
    auto totalNumInputChannels  = getTotalNumInputChannels();
    auto totalNumOutputChannels = getTotalNumOutputChannels();

    for (auto i = totalNumInputChannels; i < totalNumOutputChannels; ++i)
        buffer.clear (i, 0, buffer.getNumSamples());
    
    auto* channelDataL = buffer.getWritePointer (0);
    auto* channelDataR = buffer.getWritePointer (1);
    
    for (int i = 0; i < buffer.getNumSamples(); ++i)
    {
        for (auto membrane : membranes)
        {
            if (membrane->isExcited())
                membrane->excite();
            
            membrane->calculateFDS();
            membrane->updateStates();
            channelDataL[i] = clamp (membrane->getOutput (0.6, 0.5), -1.0, 1.0);
            channelDataR[i] = channelDataL[i];
//            std::cout << channelDataL[i] << std::endl;
        }
        
    }
}

//==============================================================================
bool MembraneAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* MembraneAudioProcessor::createEditor()
{
    return new MembraneAudioProcessorEditor (*this);
}

//==============================================================================
void MembraneAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void MembraneAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new MembraneAudioProcessor();
}

double MembraneAudioProcessor::clamp (double val, double min, double max)
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
