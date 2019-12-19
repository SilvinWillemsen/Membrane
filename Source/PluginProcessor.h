/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "Membrane.h"

//==============================================================================
/**
*/
class MembraneAudioProcessor  : public AudioProcessor
{
public:
    //==============================================================================
    MembraneAudioProcessor();
    ~MembraneAudioProcessor();

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;

   #ifndef JucePlugin_PreferredChannelConfigurations
    bool isBusesLayoutSupported (const BusesLayout& layouts) const override;
   #endif

    void processBlock (AudioBuffer<float>&, MidiBuffer&) override;

    //==============================================================================
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void changeProgramName (int index, const String& newName) override;

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;

    ////
    Membrane* getMembranePtr() { return &membrane; };
    
    double clamp (double val, double min, double max);
    
    void changeMembraneTension (int idx, double val) { membrane.setTension (val); };
    
    void changeSig0 (int idx, double val) { membrane.setSig0 (val); };
    void changeSig1 (int idx, double val) { membrane.setSig1 (val); };
    
    long getProcessorIdx() { return processorIdx; };
    
    void exciteMembrane() { membrane.setExcited(); };
    
private:
    //==============================================================================
    
    double fs = 0;
    long processorIdx = 1234;
    Membrane membrane { 400.0, 10, 0.001, 0, 15.0, 0.005, 1.0 / 44100.0, 1.0, 1.0 };
    
    AudioParameterFloat* gain;
    AudioParameterFloat* tension;
    AudioParameterFloat* sig0;
//    AudioParameterFloat* sig1;
    AudioParameterFloat* excite;
    
    bool excitedFlag = false;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MembraneAudioProcessor)
};
