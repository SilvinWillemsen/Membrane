/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"
#include "PluginProcessor.h"

//==============================================================================
/**
*/
class MembraneAudioProcessorEditor  : public AudioProcessorEditor, public Timer
{
public:
    MembraneAudioProcessorEditor (MembraneAudioProcessor&);
    ~MembraneAudioProcessorEditor();

    //==============================================================================
    void paint (Graphics&) override;
    void resized() override;
    
    void timerCallback() override;
    
private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    MembraneAudioProcessor& processor;
    std::vector<Membrane*> membranes;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MembraneAudioProcessorEditor)
};
