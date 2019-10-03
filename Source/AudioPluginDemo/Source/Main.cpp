/*
  ==============================================================================

    This file was auto-generated and contains the startup code for a PIP.

  ==============================================================================
*/

#include "../JuceLibraryCode/JuceHeader.h"
#include "AudioPluginDemo.h"

//==============================================================================
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new JuceDemoPluginAudioProcessor();
}
