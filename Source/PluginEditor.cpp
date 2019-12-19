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
    membrane = processor.getMembranePtr();

    sliders.add (new Slider());
    tensionSlider = sliders[sliders.size() - 1];
    tensionSlider->addListener (this);
    addAndMakeVisible (tensionSlider);

    sliders.add (new Slider());
    sig0Slider = sliders[sliders.size() - 1];
    sig0Slider->addListener (this);
    addAndMakeVisible (sig0Slider);

    sliders.add (new Slider());
    sig1Slider = sliders[sliders.size() - 1];
    sig1Slider->addListener (this);

    addAndMakeVisible (sig1Slider);

    sliders.add (new Slider());
    excitationWidthSlider = sliders[sliders.size() - 1];
    excitationWidthSlider->addListener (this);
    addAndMakeVisible (excitationWidthSlider);

    // specify here on which UDP port number to receive incoming OSC messages
    if (! connect (7563))                   // [3]
        showConnectionErrorMessage ("Error: could not connect to UDP port 7563.");
    
    // tell the component to listen for OSC messages matching this address:
    addListener (this, "/message"); // [4]
    
    startTimerHz (60);
    setSize (800, 600);
}

MembraneAudioProcessorEditor::~MembraneAudioProcessorEditor()
{
}

//==============================================================================
void MembraneAudioProcessorEditor::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));
    g.setColour (Colours::white);
//    exciteButton->setButtonText (String (processor.getProcessorIdx()));
    g.drawText ("Iteration " + String(iteration), getWidth() / 2.0, getHeight() / 2.0, 200, 200, Justification::centred);
}

void MembraneAudioProcessorEditor::resized()
{
    Rectangle<int> totArea = getLocalBounds();
    Rectangle<int> controlArea = totArea.removeFromBottom(getHeight() / 6.0);
    int totSlidersButtons = sliders.size() + buttons.size();
    int margin = 20;
    for (int i = 0; i < sliders.size() + buttons.size(); ++i)
    {
        if (i < sliders.size())
            sliders[i]->setBounds (controlArea.removeFromLeft (static_cast<double>(getWidth()) / static_cast<double>(totSlidersButtons)));
        else
            buttons[i - sliders.size()]->setBounds (controlArea.removeFromLeft (static_cast<int>(getWidth()) / static_cast<int>(totSlidersButtons)));
        controlArea.removeFromLeft(margin);
    }

    if (membrane != nullptr)
        membrane->setBounds (totArea);
}

void MembraneAudioProcessorEditor::timerCallback()
{
    if (updateIteration)
        iteration++;
//    double val = static_cast<double> (iteration % 255) / 127.0 - 1.0;
//    if (iteration % 100 == 0)
//        membrane->setExcited();
    if (init)
    {
        if (membrane != nullptr)
        {
            init = false;
            addAndMakeVisible (*membrane);
//            tensionSlider->setRange (1, membrane->getMaxT());
//            tensionSlider->setValue (membrane->getMaxT());
            tensionSlider->setRange (0.1, membrane->getMaxT(), 0.001);
            tensionSlider->setValue (membrane->getMaxT());
            sig0Slider->setRange (0.0, membrane->getMaxSig0());
            sig0Slider->setValue (membrane->getMaxSig0());
            sig1Slider->setRange (0.0, membrane->getMaxSig1());
            sig1Slider->setValue (membrane->getMaxSig1());
            excitationWidthSlider->setRange (2, membrane->getMaxExcitationWidth(), 1);
            excitationWidthSlider->setValue (membrane->getMaxExcitationWidth());
        }
    }
    if (membrane != nullptr)
    {
        membrane->changeDamping();
        sliders[2]->setValue (membrane->getSig1());
//        sliders[1]->setValue (membrane->getSig0());
//        sliders[0]->setValue (membrane->getTension());
    }
//    if (graphicsFlag)
        repaint();
}

void MembraneAudioProcessorEditor::sliderValueChanged(Slider* slider)
{
    if (membrane != nullptr)
    {
        if (slider == tensionSlider)
        {
//            membrane->setTensionFromSize (slider->getValue());
            membrane->setTension(slider->getValue());
//            membrane->setSig0 (membrane->getTension() / 400.0 * 15.0);
        }
        else if (slider == sig0Slider)
        {
            membrane->setSig0 (slider->getValue());
//            membrane->setTensionFromSize(slider->getValue() / 15.0);
        }
        else if (slider == sig1Slider)
        {
            membrane->setSig1 (slider->getValue());
        }
        else if (slider == excitationWidthSlider)
        {
            membrane->setExcitationWidth (slider->getValue());
        }
    }
}

void MembraneAudioProcessorEditor::buttonClicked (Button* button)
{
//    if (button == exciteButton)
//    {
//        membrane->setExcited();
//        updateIteration = !updateIteration;
//    }
//    else if (button == updateButton && membrane != nullptr)
//        membrane->setUpdate (true);
}
