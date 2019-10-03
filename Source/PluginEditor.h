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
class MembraneAudioProcessorEditor  : public AudioProcessorEditor,
                                      public Timer,
                                      public Slider::Listener,
                                      public Button::Listener
{
public:
    MembraneAudioProcessorEditor (MembraneAudioProcessor&);
    ~MembraneAudioProcessorEditor();

    //==============================================================================
    void paint (Graphics&) override;
    void resized() override;
    
    void timerCallback() override;
//
    void sliderValueChanged (Slider* slider) override;
    void buttonClicked (Button* button) override;
    
private:
    // This reference is provided as a quick way for your editor to
    // access the processor object that created it.
    MembraneAudioProcessor& processor;
    
    Membrane* membrane;
    
    OwnedArray<Slider> sliders;
    OwnedArray<TextButton> buttons;
    Slider* tensionSlider;
    Slider* sig0Slider;
    Slider* sig1Slider;
    Slider* excitationWidthSlider;

    TextButton* exciteButton;
    TextButton* updateButton;
    bool graphicsFlag = false;
    
    bool init = true;
    long iteration = 0;
    bool updateIteration = true;
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MembraneAudioProcessorEditor)
};
