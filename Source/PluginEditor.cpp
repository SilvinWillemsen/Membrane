/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin editor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"

//==============================================================================
MembraneAudioProcessorEditor::MembraneAudioProcessorEditor (MembraneAudioProcessor& p)
: AudioProcessorEditor (&p), processor (p)
{
    // Make sure that before the constructor has finished, you've set the
    // editor's size to whatever you need it to be.
    membranes.push_back (processor.getMembranePtr (0));
    addAndMakeVisible (membranes[0]);
    startTimerHz (60);
    setSize (400, 300);
}

MembraneAudioProcessorEditor::~MembraneAudioProcessorEditor()
{
}

//==============================================================================
void MembraneAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));
}

void MembraneAudioProcessorEditor::resized()
{
    membranes[0]->setBounds(getLocalBounds());
}

void MembraneAudioProcessorEditor::timerCallback()
{
    repaint();
}
