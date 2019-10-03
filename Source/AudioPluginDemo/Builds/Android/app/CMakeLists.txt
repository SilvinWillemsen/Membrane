# Automatically generated makefile, created by the Projucer
# Don't edit this file! Your changes will be overwritten when you re-save the Projucer project!

cmake_minimum_required(VERSION 3.4.1)

SET(BINARY_NAME "juce_jni")

add_library("cpufeatures" STATIC "${ANDROID_NDK}/sources/android/cpufeatures/cpu-features.c")
set_source_files_properties("${ANDROID_NDK}/sources/android/cpufeatures/cpu-features.c" PROPERTIES COMPILE_FLAGS "-Wno-sign-conversion -Wno-gnu-statement-expression")

add_definitions("-DJUCE_ANDROID=1" "-DJUCE_ANDROID_API_VERSION=10" "-DJUCE_ANDROID_ACTIVITY_CLASSNAME=com_juce_audioplugindemo_AudioPluginDemo" "-DJUCE_ANDROID_ACTIVITY_CLASSPATH=\"com/juce/audioplugindemo/AudioPluginDemo\"" "-DJUCE_ANDROID_SHARING_CONTENT_PROVIDER_CLASSNAME=com_juce_audioplugindemo_SharingContentProvider" "-DJUCE_ANDROID_SHARING_CONTENT_PROVIDER_CLASSPATH=\"com/juce/audioplugindemo/SharingContentProvider\"" "-DJUCE_PUSH_NOTIFICATIONS=1" "-DPIP_JUCE_EXAMPLES_DIRECTORY=L1VzZXJzL1NpbHZpblcvcmVwb3NpdG9yaWVzL25ld0pVQ0UvSlVDRS9leGFtcGxlcw==" "-DJUCER_ANDROIDSTUDIO_7F0E4A25=1" "-DJUCE_APP_VERSION=1.0.0" "-DJUCE_APP_VERSION_HEX=0x10000")

include_directories( AFTER
    "../../../JuceLibraryCode"
    "/Users/SilvinW/repositories/newJUCE/JUCE/modules"
    "${ANDROID_NDK}/sources/android/cpufeatures"
)

enable_language(ASM)

IF(JUCE_BUILD_CONFIGURATION MATCHES "DEBUG")
    add_definitions("-DDEBUG=1" "-D_DEBUG=1")
ELSEIF(JUCE_BUILD_CONFIGURATION MATCHES "RELEASE")
    add_definitions("-DNDEBUG=1")
ELSE(JUCE_BUILD_CONFIGURATION MATCHES "DEBUG")
    MESSAGE( FATAL_ERROR "No matching build-configuration found." )
ENDIF(JUCE_BUILD_CONFIGURATION MATCHES "DEBUG")

add_library( ${BINARY_NAME}

    SHARED

    "../../../Source/Main.cpp"
    "../../../Source/AudioPluginDemo.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/audio_play_head/juce_AudioPlayHead.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioChannelSet.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioChannelSet.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioDataConverters.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioDataConverters.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioProcessLoadMeasurer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioProcessLoadMeasurer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioSampleBuffer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_FloatVectorOperations.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_FloatVectorOperations.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_CatmullRomInterpolator.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_CatmullRomInterpolator.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_Decibels.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_IIRFilter.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_IIRFilter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_LagrangeInterpolator.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_LagrangeInterpolator.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_LinearSmoothedValue.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_Reverb.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiBuffer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiBuffer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiFile.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiFile.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiKeyboardState.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiKeyboardState.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiMessage.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiMessage.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiMessageSequence.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiMessageSequence.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiRPN.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiRPN.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEInstrument.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEInstrument.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEMessages.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEMessages.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPENote.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPENote.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiser.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiserBase.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiserBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiserVoice.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiserVoice.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEUtils.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEUtils.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEValue.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEValue.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEZoneLayout.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEZoneLayout.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/native/juce_mac_CoreAudioLayouts.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_AudioSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_BufferingAudioSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_BufferingAudioSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ChannelRemappingAudioSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ChannelRemappingAudioSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_IIRFilterAudioSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_IIRFilterAudioSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_MemoryAudioSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_MemoryAudioSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_MixerAudioSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_MixerAudioSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_PositionableAudioSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ResamplingAudioSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ResamplingAudioSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ReverbAudioSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ReverbAudioSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ToneGeneratorAudioSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ToneGeneratorAudioSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/synthesisers/juce_Synthesiser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/synthesisers/juce_Synthesiser.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/juce_audio_basics.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/juce_audio_basics.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_basics/juce_audio_basics.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioDeviceManager.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioDeviceManager.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioIODevice.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioIODevice.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioIODeviceType.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioIODeviceType.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_SystemAudioVolume.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/midi_io/juce_MidiInput.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/midi_io/juce_MidiMessageCollector.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/midi_io/juce_MidiMessageCollector.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/midi_io/juce_MidiOutput.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/midi_io/juce_MidiOutput.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_android_Audio.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_android_Midi.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_android_Oboe.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_android_OpenSL.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_ios_Audio.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_ios_Audio.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_linux_ALSA.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_linux_Bela.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_linux_JackAudio.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_linux_Midi.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_mac_CoreAudio.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_mac_CoreMidi.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_MidiDataConcatenator.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_win32_ASIO.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_win32_DirectSound.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_win32_Midi.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_win32_WASAPI.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/sources/juce_AudioSourcePlayer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/sources/juce_AudioSourcePlayer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/sources/juce_AudioTransportSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/sources/juce_AudioTransportSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/juce_audio_devices.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/juce_audio_devices.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_devices/juce_audio_devices.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/all.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/bitmath.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/bitreader.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/bitwriter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/cpu.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/crc.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/fixed.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/float.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/format.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/lpc.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/md5.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/memory.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/metadata.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/stream_encoder.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/stream_encoder_framing.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/window.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/protected/all.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/protected/stream_decoder.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/protected/stream_encoder.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/bitmath.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/bitreader.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/bitwriter.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/cpu.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/crc.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/fixed.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/float.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/format.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/lpc_flac.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/md5.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/memory.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/stream_decoder.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/stream_encoder.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/stream_encoder_framing.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/window_flac.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/all.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/alloc.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/assert.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/callback.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/compat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/endswap.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/export.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/Flac Licence.txt"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/format.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/metadata.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/ordinals.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/stream_decoder.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/stream_encoder.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/win_utf8_io.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/books/coupled/res_books_51.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/books/coupled/res_books_stereo.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/books/floor/floor_books.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/books/uncoupled/res_books_uncoupled.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/floor_all.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/psych_8.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/psych_11.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/psych_16.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/psych_44.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/residue_8.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/residue_16.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/residue_44.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/residue_44p51.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/residue_44u.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_8.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_11.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_16.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_22.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_32.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_44.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_44p51.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_44u.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_X.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/analysis.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/backends.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/bitrate.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/bitrate.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/block.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/codebook.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/codebook.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/codec_internal.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/envelope.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/envelope.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/floor0.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/floor1.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/highlevel.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/info.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lookup.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lookup.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lookup_data.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lpc.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lpc.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lsp.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lsp.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/mapping0.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/masking.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/mdct.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/mdct.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/misc.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/os.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/psy.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/psy.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/registry.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/registry.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/res0.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/scales.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/sharedbook.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/smallft.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/smallft.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/synthesis.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/vorbisenc.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/vorbisfile.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/window.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/window.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/bitwise.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/codec.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/config_types.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/framing.c"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/Ogg Vorbis Licence.txt"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/ogg.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/os_types.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/vorbisenc.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/vorbisfile.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_AiffAudioFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_AiffAudioFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_CoreAudioFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_CoreAudioFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_FlacAudioFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_FlacAudioFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_LAMEEncoderAudioFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_LAMEEncoderAudioFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_MP3AudioFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_MP3AudioFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_OggVorbisAudioFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_OggVorbisAudioFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_WavAudioFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_WavAudioFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_WindowsMediaAudioFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_WindowsMediaAudioFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatManager.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatManager.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatReader.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatReader.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatReaderSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatReaderSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatWriter.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatWriter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioSubsectionReader.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioSubsectionReader.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_BufferingAudioFormatReader.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_BufferingAudioFormatReader.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_MemoryMappedAudioFormatReader.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/sampler/juce_Sampler.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/sampler/juce_Sampler.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/juce_audio_formats.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/juce_audio_formats.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_formats/juce_audio_formats.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AAX/juce_AAX_Modifier_Injector.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AAX/juce_AAX_Wrapper.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUBase.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUBaseHelper.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUBuffer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUBuffer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUCarbonViewBase.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUCarbonViewBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUCarbonViewControl.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUCarbonViewControl.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUCarbonViewDispatch.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUDispatch.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUDispatch.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUInputElement.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUInputElement.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUInputFormatConverter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUMIDIBase.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUMIDIBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUOutputBase.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUOutputBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUOutputElement.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUOutputElement.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUPlugInDispatch.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUPlugInDispatch.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUScopeElement.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUScopeElement.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUSilentTimeout.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUTimestampGenerator.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUViewLocalizedStringKeys.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAtomic.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAtomicStack.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAudioChannelLayout.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAudioChannelLayout.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAUParameter.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAUParameter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAutoDisposer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CADebugMacros.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CADebugPrintf.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAException.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAHostTimeBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAMath.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAMutex.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAMutex.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CarbonEventHandler.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CarbonEventHandler.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAReferenceCounted.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAStreamBasicDescription.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAStreamBasicDescription.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAThreadSafeList.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAVectorUnit.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAVectorUnit.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAVectorUnitTypes.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAXException.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/ComponentBase.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/ComponentBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/MusicDeviceBase.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/MusicDeviceBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/juce_AU_Wrapper.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/juce_AUv3_Wrapper.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_DigiCode_Header.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_DigiCode1.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_DigiCode2.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_DigiCode3.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_MacUtilities.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_WinUtilities.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_Wrapper.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/Standalone/juce_StandaloneFilterApp.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/Standalone/juce_StandaloneFilterWindow.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/Unity/juce_Unity_Wrapper.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/Unity/juce_UnityPluginInterface.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_CarbonVisibility.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_CheckSettingMacros.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_FakeMouseMoveGenerator.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_IncludeModuleHeaders.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_IncludeSystemHeaders.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_PluginHostType.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_PluginUtilities.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_WindowsHooks.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/VST/juce_VST_Wrapper.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/VST/juce_VST_Wrapper.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/VST/juce_VSTCallbackHandler.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/VST3/juce_VST3_Wrapper.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_AAX.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_AAX.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_AU_1.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_AU_2.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_AUv3.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_1.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_2.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_3.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_4.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_utils.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_utils.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_Standalone.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_Unity.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_utils.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_VST_utils.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_VST2.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_VST3.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format/juce_AudioPluginFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format/juce_AudioPluginFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format/juce_AudioPluginFormatManager.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format/juce_AudioPluginFormatManager.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/baseiids.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/classfactoryhelpers.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fbuffer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fbuffer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fdebug.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fdebug.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fobject.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fobject.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fstreamer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fstreamer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fstring.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fstring.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/updatehandler.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/updatehandler.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/thread/include/flock.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/thread/source/flock.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/LICENSE.txt"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/conststringtable.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/conststringtable.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/coreiids.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/falignpop.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/falignpush.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/fplatform.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/fstrdefs.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ftypes.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/funknown.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/funknown.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/futils.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/fvariant.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ibstream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/icloneable.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ipersistent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ipluginbase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/istringresult.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/iupdatehandler.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/smartpointer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ustring.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ustring.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/gui/iplugview.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/gui/iplugviewcontentscalesupport.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstattributes.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstaudioprocessor.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstautomationstate.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstchannelcontextinfo.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstcomponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstcontextmenu.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivsteditcontroller.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstevents.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivsthostapplication.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstinterappaudio.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstmessage.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstmidicontrollers.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstnoteexpression.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstparameterchanges.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstplugview.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstprefetchablesupport.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstprocesscontext.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstrepresentation.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstunits.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/vstpshpack4.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/vstspeaker.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/vsttypes.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/LICENSE.txt"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/common/memorystream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/common/memorystream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/common/pluginview.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/common/pluginview.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/hosting/hostclasses.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/hosting/hostclasses.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstbus.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstbus.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstcomponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstcomponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstcomponentbase.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstcomponentbase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vsteditcontroller.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vsteditcontroller.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstinitiids.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstparameters.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstparameters.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstpresetfile.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstpresetfile.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/LICENSE.txt"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/LICENSE.txt"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/README.md"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_AU_Shared.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_AudioUnitPluginFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_AudioUnitPluginFormat.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_LADSPAPluginFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_LADSPAPluginFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_LegacyAudioParameter.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VST3Common.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VST3Headers.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VST3PluginFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VST3PluginFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VSTCommon.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VSTMidiEventList.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VSTPluginFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VSTPluginFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioPluginInstance.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioPluginInstance.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessor.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessor.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorEditor.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorEditor.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorGraph.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorGraph.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorListener.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorParameter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorParameterGroup.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorParameterGroup.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_GenericAudioProcessorEditor.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_GenericAudioProcessorEditor.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_PluginDescription.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_PluginDescription.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_KnownPluginList.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_KnownPluginList.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_PluginDirectoryScanner.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_PluginDirectoryScanner.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_PluginListComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_PluginListComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioParameterBool.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioParameterChoice.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioParameterFloat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioParameterInt.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioProcessorParameters.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioProcessorParameterWithID.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioProcessorValueTreeState.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioProcessorValueTreeState.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_RangedAudioParameter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/juce_audio_processors.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/juce_audio_processors.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_processors/juce_audio_processors.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/audio_cd/juce_AudioCDBurner.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/audio_cd/juce_AudioCDReader.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/audio_cd/juce_AudioCDReader.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioAppComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioAppComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioDeviceSelectorComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioDeviceSelectorComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioThumbnail.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioThumbnail.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioThumbnailBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioThumbnailCache.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioThumbnailCache.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioVisualiserComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioVisualiserComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_BluetoothMidiDevicePairingDialogue.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_MidiKeyboardComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_MidiKeyboardComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_android_BluetoothMidiDevicePairingDialogue.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_ios_BluetoothMidiDevicePairingDialogue.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_linux_AudioCDReader.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_linux_BluetoothMidiDevicePairingDialogue.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_mac_AudioCDBurner.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_mac_AudioCDReader.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_mac_BluetoothMidiDevicePairingDialogue.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_win_BluetoothMidiDevicePairingDialogue.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_win32_AudioCDBurner.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_win32_AudioCDReader.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/players/juce_AudioProcessorPlayer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/players/juce_AudioProcessorPlayer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/players/juce_SoundPlayer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/players/juce_SoundPlayer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/juce_audio_utils.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/juce_audio_utils.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_audio_utils/juce_audio_utils.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_AbstractFifo.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_AbstractFifo.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_Array.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ArrayAllocationBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ArrayBase.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ArrayBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_DynamicObject.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_DynamicObject.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ElementComparator.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_HashMap.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_HashMap_test.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_LinkedListPointer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ListenerList.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_NamedValueSet.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_NamedValueSet.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_OwnedArray.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_OwnedArray.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_PropertySet.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_PropertySet.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ReferenceCountedArray.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ReferenceCountedArray.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ScopedValueSetter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_SortedSet.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_SparseSet.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_SparseSet.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_Variant.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_Variant.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_DirectoryIterator.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_DirectoryIterator.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_File.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_File.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileFilter.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileFilter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileInputStream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileInputStream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileOutputStream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileOutputStream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileSearchPath.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileSearchPath.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_MemoryMappedFile.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_TemporaryFile.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_TemporaryFile.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_WildcardFileFilter.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_WildcardFileFilter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/javascript/juce_Javascript.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/javascript/juce_Javascript.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/javascript/juce_JSON.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/javascript/juce_JSON.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/logging/juce_FileLogger.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/logging/juce_FileLogger.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/logging/juce_Logger.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/logging/juce_Logger.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_BigInteger.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_BigInteger.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_Expression.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_Expression.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_MathsFunctions.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_NormalisableRange.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_Random.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_Random.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_Range.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_StatisticsAccumulator.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_Atomic.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_ByteOrder.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_ContainerDeletePolicy.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_HeapBlock.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_HeavyweightLeakedObjectDetector.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_LeakedObjectDetector.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_Memory.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_MemoryBlock.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_MemoryBlock.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_OptionalScopedPointer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_ReferenceCountedObject.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_ScopedPointer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_SharedResourcePointer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_Singleton.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_WeakReference.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_ConsoleApplication.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_ConsoleApplication.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_Result.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_Result.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_RuntimePermissions.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_RuntimePermissions.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_StdFunctionCompat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_StdFunctionCompat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_Uuid.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_Uuid.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_WindowsRegistry.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_Files.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_JNIHelpers.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_Misc.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_Network.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_RuntimePermissions.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_SystemStats.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_Threads.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_BasicNativeHeaders.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_curl_Network.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_linux_CommonFile.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_linux_Files.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_linux_Network.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_linux_SystemStats.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_linux_Threads.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_ClangBugWorkaround.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_Files.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_Network.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_Strings.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_SystemStats.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_Threads.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_osx_ObjCHelpers.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_posix_IPAddress.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_posix_NamedPipe.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_posix_SharedCode.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_ComSmartPtr.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_Files.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_Network.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_Registry.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_SystemStats.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_Threads.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_IPAddress.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_IPAddress.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_MACAddress.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_MACAddress.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_NamedPipe.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_NamedPipe.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_Socket.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_Socket.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_URL.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_URL.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_WebInputStream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_WebInputStream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_BufferedInputStream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_BufferedInputStream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_FileInputSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_FileInputSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_InputSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_InputStream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_InputStream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_MemoryInputStream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_MemoryInputStream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_MemoryOutputStream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_MemoryOutputStream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_OutputStream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_OutputStream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_SubregionStream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_SubregionStream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_URLInputSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_URLInputSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_CompilerSupport.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_PlatformDefs.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_StandardHeader.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_SystemStats.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_SystemStats.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_TargetPlatform.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_Base64.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_Base64.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharacterFunctions.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharacterFunctions.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharPointer_ASCII.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharPointer_UTF8.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharPointer_UTF16.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharPointer_UTF32.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_Identifier.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_Identifier.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_LocalisedStrings.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_LocalisedStrings.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_NewLine.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_String.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_String.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringArray.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringArray.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringPairArray.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringPairArray.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringPool.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringPool.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringRef.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_TextDiff.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_TextDiff.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ChildProcess.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ChildProcess.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_CriticalSection.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_DynamicLibrary.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_HighResolutionTimer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_HighResolutionTimer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_InterProcessLock.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_Process.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ReadWriteLock.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ReadWriteLock.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ScopedLock.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ScopedReadLock.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ScopedWriteLock.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_SpinLock.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_Thread.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_Thread.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ThreadLocalValue.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ThreadPool.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ThreadPool.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_TimeSliceThread.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_TimeSliceThread.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_WaitableEvent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_PerformanceCounter.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_PerformanceCounter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_RelativeTime.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_RelativeTime.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_Time.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_Time.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/unit_tests/juce_UnitTest.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/unit_tests/juce_UnitTest.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/xml/juce_XmlDocument.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/xml/juce_XmlDocument.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/xml/juce_XmlElement.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/xml/juce_XmlElement.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/adler32.c"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/compress.c"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/crc32.c"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/crc32.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/deflate.c"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/deflate.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/infback.c"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inffast.c"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inffast.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inffixed.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inflate.c"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inflate.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inftrees.c"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inftrees.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/trees.c"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/trees.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/uncompr.c"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/zconf.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/zconf.in.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/zlib.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/zutil.c"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/zutil.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_GZIPCompressorOutputStream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_GZIPCompressorOutputStream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_GZIPDecompressorInputStream.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_GZIPDecompressorInputStream.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_ZipFile.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_ZipFile.h"
    "../../../../../../newJUCE/JUCE/modules/juce_core/juce_core.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_core/juce_core.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_core/juce_core.h"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/app_properties/juce_ApplicationProperties.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/app_properties/juce_ApplicationProperties.h"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/app_properties/juce_PropertiesFile.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/app_properties/juce_PropertiesFile.h"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/undomanager/juce_UndoableAction.h"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/undomanager/juce_UndoManager.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/undomanager/juce_UndoManager.h"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_CachedValue.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_CachedValue.h"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_Value.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_Value.h"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_ValueTree.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_ValueTree.h"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_ValueTreeSynchroniser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_ValueTreeSynchroniser.h"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_ValueWithDefault.h"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/juce_data_structures.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/juce_data_structures.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_data_structures/juce_data_structures.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ActionBroadcaster.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ActionBroadcaster.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ActionListener.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_AsyncUpdater.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_AsyncUpdater.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ChangeBroadcaster.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ChangeBroadcaster.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ChangeListener.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_ConnectedChildProcess.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_ConnectedChildProcess.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_InterprocessConnection.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_InterprocessConnection.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_InterprocessConnectionServer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_InterprocessConnectionServer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_NetworkServiceDiscovery.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_NetworkServiceDiscovery.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_ApplicationBase.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_ApplicationBase.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_CallbackMessage.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_DeletedAtShutdown.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_DeletedAtShutdown.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_Initialisation.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_Message.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_MessageListener.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_MessageListener.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_MessageManager.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_MessageManager.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_MountedVolumeListChangeDetector.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_NotificationType.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_android_Messaging.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_ios_MessageManager.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_linux_EventLoop.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_linux_Messaging.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_mac_MessageManager.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_osx_MessageQueue.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_win32_HiddenMessageWindow.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_win32_Messaging.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_win32_WinRTWrapper.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_win32_WinRTWrapper.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/timers/juce_MultiTimer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/timers/juce_MultiTimer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/timers/juce_Timer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/timers/juce_Timer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_events/juce_events.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_events/juce_events.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_events/juce_events.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_Colour.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_Colour.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_ColourGradient.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_ColourGradient.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_Colours.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_Colours.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_FillType.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_FillType.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_PixelFormats.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_GraphicsContext.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_GraphicsContext.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_LowLevelGraphicsContext.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_LowLevelGraphicsPostScriptRenderer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_LowLevelGraphicsPostScriptRenderer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_LowLevelGraphicsSoftwareRenderer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_LowLevelGraphicsSoftwareRenderer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/effects/juce_DropShadowEffect.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/effects/juce_DropShadowEffect.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/effects/juce_GlowEffect.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/effects/juce_GlowEffect.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/effects/juce_ImageEffectFilter.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_AttributedString.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_AttributedString.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_CustomTypeface.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_CustomTypeface.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_Font.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_Font.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_GlyphArrangement.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_GlyphArrangement.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_TextLayout.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_TextLayout.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_Typeface.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_Typeface.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_AffineTransform.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_AffineTransform.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_BorderSize.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_EdgeTable.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_EdgeTable.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Line.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Parallelogram.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Path.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Path.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_PathIterator.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_PathIterator.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_PathStrokeType.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_PathStrokeType.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Point.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Rectangle.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_RectangleList.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/cderror.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/changes to libjpeg for JUCE.txt"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcapimin.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcapistd.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jccoefct.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jccolor.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcdctmgr.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jchuff.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jchuff.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcinit.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcmainct.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcmarker.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcmaster.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcomapi.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jconfig.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcparam.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcphuff.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcprepct.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcsample.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jctrans.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdapimin.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdapistd.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdatasrc.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdcoefct.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdcolor.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdct.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jddctmgr.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdhuff.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdhuff.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdinput.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdmainct.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdmarker.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdmaster.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdmerge.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdphuff.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdpostct.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdsample.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdtrans.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jerror.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jerror.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jfdctflt.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jfdctfst.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jfdctint.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jidctflt.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jidctfst.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jidctint.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jidctred.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jinclude.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jmemmgr.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jmemnobs.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jmemsys.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jmorecfg.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jpegint.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jpeglib.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jquant1.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jquant2.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jutils.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jversion.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/transupp.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/transupp.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/libpng_readme.txt"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/png.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/png.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngconf.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngerror.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngget.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pnginfo.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngmem.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngpread.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngpriv.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngread.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngrio.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngrtran.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngrutil.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngset.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngstruct.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngtrans.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngwio.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngwrite.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngwtran.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngwutil.c"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/juce_GIFLoader.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/juce_JPEGLoader.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/juce_PNGLoader.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_Image.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_Image.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageCache.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageCache.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageConvolutionKernel.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageConvolutionKernel.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageFileFormat.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageFileFormat.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_android_Fonts.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_android_GraphicsContext.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_android_IconHelpers.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_freetype_Fonts.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_linux_Fonts.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_linux_IconHelpers.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_mac_CoreGraphicsContext.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_mac_CoreGraphicsContext.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_mac_CoreGraphicsHelpers.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_mac_Fonts.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_mac_IconHelpers.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_RenderingHelpers.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_Direct2DGraphicsContext.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_Direct2DGraphicsContext.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_DirectWriteTypeface.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_DirectWriteTypeLayout.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_Fonts.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_IconHelpers.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/placement/juce_Justification.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/placement/juce_RectanglePlacement.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/placement/juce_RectanglePlacement.h"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/juce_graphics.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/juce_graphics.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_graphics/juce_graphics.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/application/juce_Application.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/application/juce_Application.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ArrowButton.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ArrowButton.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_Button.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_Button.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_DrawableButton.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_DrawableButton.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_HyperlinkButton.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_HyperlinkButton.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ImageButton.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ImageButton.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ShapeButton.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ShapeButton.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_TextButton.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_TextButton.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ToggleButton.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ToggleButton.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ToolbarButton.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ToolbarButton.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandID.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandInfo.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandInfo.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandManager.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandManager.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandTarget.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandTarget.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_KeyPressMappingSet.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_KeyPressMappingSet.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_CachedComponentImage.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_Component.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_Component.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_ComponentListener.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_ComponentListener.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_ModalComponentManager.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_ModalComponentManager.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/desktop/juce_Desktop.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/desktop/juce_Desktop.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/desktop/juce_Displays.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/desktop/juce_Displays.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_Drawable.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_Drawable.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableComposite.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableComposite.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableImage.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableImage.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawablePath.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawablePath.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableRectangle.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableRectangle.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableShape.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableShape.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableText.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableText.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_SVGParser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_ContentSharer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_ContentSharer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_DirectoryContentsDisplayComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_DirectoryContentsDisplayComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_DirectoryContentsList.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_DirectoryContentsList.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileBrowserComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileBrowserComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileBrowserListener.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileChooser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileChooser.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileChooserDialogBox.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileChooserDialogBox.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileListComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileListComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FilenameComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FilenameComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FilePreviewComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileSearchPathListComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileSearchPathListComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileTreeComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileTreeComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_ImagePreviewComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_ImagePreviewComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_CaretComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_CaretComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyboardFocusTraverser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyboardFocusTraverser.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyListener.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyListener.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyPress.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyPress.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_ModifierKeys.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_ModifierKeys.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_SystemClipboard.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_TextEditorKeyMapper.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_TextInputTarget.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_AnimatedPosition.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_AnimatedPositionBehaviours.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentAnimator.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentAnimator.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentBoundsConstrainer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentBoundsConstrainer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentBuilder.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentBuilder.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentMovementWatcher.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentMovementWatcher.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ConcertinaPanel.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ConcertinaPanel.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_FlexBox.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_FlexBox.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_FlexItem.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_Grid.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_Grid.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_GridItem.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_GridItem.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_GridUnitTests.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_GroupComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_GroupComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_MultiDocumentPanel.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_MultiDocumentPanel.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableBorderComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableBorderComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableCornerComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableCornerComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableEdgeComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableEdgeComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ScrollBar.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ScrollBar.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_SidePanel.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_SidePanel.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableLayoutManager.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableLayoutManager.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableLayoutResizerBar.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableLayoutResizerBar.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableObjectResizer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableObjectResizer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_TabbedButtonBar.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_TabbedButtonBar.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_TabbedComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_TabbedComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_Viewport.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_Viewport.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V1.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V1.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V2.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V2.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V3.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V3.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V4.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V4.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_BurgerMenuComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_BurgerMenuComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_MenuBarComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_MenuBarComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_MenuBarModel.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_MenuBarModel.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_PopupMenu.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_PopupMenu.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_BubbleComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_BubbleComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_DropShadower.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_DropShadower.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_JUCESplashScreen.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_JUCESplashScreen.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_ComponentDragger.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_ComponentDragger.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_DragAndDropContainer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_DragAndDropContainer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_DragAndDropTarget.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_FileDragAndDropTarget.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_LassoComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseCursor.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseCursor.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseEvent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseEvent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseInactivityDetector.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseInactivityDetector.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseInputSource.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseInputSource.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseListener.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseListener.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_SelectedItemSet.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_TextDragAndDropTarget.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_TooltipClient.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_android_ContentSharer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_android_FileChooser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_android_Windowing.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_common_MimeTypes.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_ios_ContentSharer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_ios_FileChooser.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_ios_UIViewComponentPeer.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_ios_Windowing.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_linux_FileChooser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_linux_X11.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_linux_X11.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_linux_X11_Clipboard.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_linux_X11_Windowing.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_mac_FileChooser.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_mac_MainMenu.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_mac_MouseCursor.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_mac_NSViewComponentPeer.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_mac_Windowing.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_MultiTouchMapper.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_win32_DragAndDrop.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_win32_FileChooser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_win32_Windowing.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_MarkerList.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_MarkerList.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeCoordinate.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeCoordinate.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeCoordinatePositioner.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeCoordinatePositioner.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeParallelogram.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeParallelogram.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativePoint.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativePoint.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativePointPath.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativePointPath.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeRectangle.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeRectangle.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_BooleanPropertyComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_BooleanPropertyComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_ButtonPropertyComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_ButtonPropertyComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_ChoicePropertyComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_ChoicePropertyComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_MultiChoicePropertyComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_MultiChoicePropertyComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_PropertyComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_PropertyComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_PropertyPanel.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_PropertyPanel.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_SliderPropertyComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_SliderPropertyComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_TextPropertyComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_TextPropertyComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ComboBox.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ComboBox.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ImageComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ImageComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Label.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Label.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ListBox.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ListBox.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ProgressBar.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ProgressBar.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Slider.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Slider.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TableHeaderComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TableHeaderComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TableListBox.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TableListBox.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TextEditor.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TextEditor.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Toolbar.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Toolbar.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ToolbarItemComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ToolbarItemComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ToolbarItemFactory.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ToolbarItemPalette.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ToolbarItemPalette.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TreeView.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TreeView.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_AlertWindow.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_AlertWindow.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_CallOutBox.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_CallOutBox.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ComponentPeer.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ComponentPeer.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_DialogWindow.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_DialogWindow.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_DocumentWindow.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_DocumentWindow.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_NativeMessageBox.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ResizableWindow.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ResizableWindow.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ThreadWithProgressWindow.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ThreadWithProgressWindow.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_TooltipWindow.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_TooltipWindow.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_TopLevelWindow.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_TopLevelWindow.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/juce_gui_basics.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/juce_gui_basics.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_basics/juce_gui_basics.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CodeDocument.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CodeDocument.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CodeEditorComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CodeEditorComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CodeTokeniser.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CPlusPlusCodeTokeniser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CPlusPlusCodeTokeniser.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CPlusPlusCodeTokeniserFunctions.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_LuaCodeTokeniser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_LuaCodeTokeniser.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_XMLCodeTokeniser.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_XMLCodeTokeniser.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/documents/juce_FileBasedDocument.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/documents/juce_FileBasedDocument.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_ActiveXControlComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_AndroidViewComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_NSViewComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_ScopedDPIAwarenessDisabler.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_UIViewComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_XEmbedComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_AnimatedAppComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_AnimatedAppComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_AppleRemote.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_BubbleMessageComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_BubbleMessageComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_ColourSelector.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_ColourSelector.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_KeyMappingEditorComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_KeyMappingEditorComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_LiveConstantEditor.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_LiveConstantEditor.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_PreferencesPanel.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_PreferencesPanel.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_PushNotifications.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_PushNotifications.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_RecentlyOpenedFilesList.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_RecentlyOpenedFilesList.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_SplashScreen.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_SplashScreen.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_SystemTrayIconComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_SystemTrayIconComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_WebBrowserComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_android_PushNotifications.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_android_WebBrowserComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_AndroidViewComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_ios_PushNotifications.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_ios_UIViewComponent.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_linux_X11_SystemTrayIcon.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_linux_X11_WebBrowserComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_linux_XEmbedComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_AppleRemote.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_CarbonViewWrapperComponent.h"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_NSViewComponent.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_PushNotifications.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_SystemTrayIcon.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_WebBrowserComponent.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_win32_ActiveXComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_win32_SystemTrayIcon.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_win32_WebBrowserComponent.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/juce_gui_extra.cpp"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/juce_gui_extra.mm"
    "../../../../../../newJUCE/JUCE/modules/juce_gui_extra/juce_gui_extra.h"
    "../../../JuceLibraryCode/AppConfig.h"
    "../../../JuceLibraryCode/include_juce_audio_basics.cpp"
    "../../../JuceLibraryCode/include_juce_audio_devices.cpp"
    "../../../JuceLibraryCode/include_juce_audio_formats.cpp"
    "../../../JuceLibraryCode/include_juce_audio_plugin_client_Standalone.cpp"
    "../../../JuceLibraryCode/include_juce_audio_plugin_client_utils.cpp"
    "../../../JuceLibraryCode/include_juce_audio_processors.cpp"
    "../../../JuceLibraryCode/include_juce_audio_utils.cpp"
    "../../../JuceLibraryCode/include_juce_core.cpp"
    "../../../JuceLibraryCode/include_juce_data_structures.cpp"
    "../../../JuceLibraryCode/include_juce_events.cpp"
    "../../../JuceLibraryCode/include_juce_graphics.cpp"
    "../../../JuceLibraryCode/include_juce_gui_basics.cpp"
    "../../../JuceLibraryCode/include_juce_gui_extra.cpp"
    "../../../JuceLibraryCode/JuceHeader.h"
)

set_source_files_properties("../../../Source/AudioPluginDemo.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/audio_play_head/juce_AudioPlayHead.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioChannelSet.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioChannelSet.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioDataConverters.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioDataConverters.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioProcessLoadMeasurer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioProcessLoadMeasurer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_AudioSampleBuffer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_FloatVectorOperations.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/buffers/juce_FloatVectorOperations.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_CatmullRomInterpolator.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_CatmullRomInterpolator.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_Decibels.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_IIRFilter.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_IIRFilter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_LagrangeInterpolator.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_LagrangeInterpolator.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_LinearSmoothedValue.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/effects/juce_Reverb.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiBuffer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiBuffer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiFile.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiFile.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiKeyboardState.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiKeyboardState.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiMessage.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiMessage.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiMessageSequence.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiMessageSequence.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiRPN.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/midi/juce_MidiRPN.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEInstrument.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEInstrument.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEMessages.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEMessages.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPENote.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPENote.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiser.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiserBase.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiserBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiserVoice.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPESynthesiserVoice.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEUtils.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEUtils.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEValue.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEValue.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEZoneLayout.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/mpe/juce_MPEZoneLayout.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/native/juce_mac_CoreAudioLayouts.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_AudioSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_BufferingAudioSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_BufferingAudioSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ChannelRemappingAudioSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ChannelRemappingAudioSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_IIRFilterAudioSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_IIRFilterAudioSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_MemoryAudioSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_MemoryAudioSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_MixerAudioSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_MixerAudioSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_PositionableAudioSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ResamplingAudioSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ResamplingAudioSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ReverbAudioSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ReverbAudioSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ToneGeneratorAudioSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/sources/juce_ToneGeneratorAudioSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/synthesisers/juce_Synthesiser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/synthesisers/juce_Synthesiser.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/juce_audio_basics.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/juce_audio_basics.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_basics/juce_audio_basics.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioDeviceManager.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioDeviceManager.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioIODevice.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioIODevice.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioIODeviceType.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_AudioIODeviceType.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/audio_io/juce_SystemAudioVolume.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/midi_io/juce_MidiInput.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/midi_io/juce_MidiMessageCollector.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/midi_io/juce_MidiMessageCollector.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/midi_io/juce_MidiOutput.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/midi_io/juce_MidiOutput.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_android_Audio.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_android_Midi.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_android_Oboe.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_android_OpenSL.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_ios_Audio.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_ios_Audio.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_linux_ALSA.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_linux_Bela.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_linux_JackAudio.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_linux_Midi.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_mac_CoreAudio.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_mac_CoreMidi.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_MidiDataConcatenator.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_win32_ASIO.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_win32_DirectSound.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_win32_Midi.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/native/juce_win32_WASAPI.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/sources/juce_AudioSourcePlayer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/sources/juce_AudioSourcePlayer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/sources/juce_AudioTransportSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/sources/juce_AudioTransportSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/juce_audio_devices.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/juce_audio_devices.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_devices/juce_audio_devices.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/all.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/bitmath.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/bitreader.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/bitwriter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/cpu.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/crc.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/fixed.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/float.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/format.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/lpc.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/md5.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/memory.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/metadata.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/stream_encoder.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/stream_encoder_framing.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/private/window.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/protected/all.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/protected/stream_decoder.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/include/protected/stream_encoder.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/bitmath.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/bitreader.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/bitwriter.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/cpu.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/crc.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/fixed.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/float.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/format.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/lpc_flac.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/md5.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/memory.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/stream_decoder.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/stream_encoder.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/stream_encoder_framing.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/libFLAC/window_flac.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/all.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/alloc.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/assert.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/callback.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/compat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/endswap.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/export.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/Flac Licence.txt" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/format.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/metadata.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/ordinals.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/stream_decoder.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/stream_encoder.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/flac/win_utf8_io.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/books/coupled/res_books_51.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/books/coupled/res_books_stereo.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/books/floor/floor_books.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/books/uncoupled/res_books_uncoupled.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/floor_all.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/psych_8.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/psych_11.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/psych_16.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/psych_44.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/residue_8.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/residue_16.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/residue_44.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/residue_44p51.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/residue_44u.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_8.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_11.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_16.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_22.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_32.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_44.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_44p51.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_44u.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/modes/setup_X.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/analysis.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/backends.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/bitrate.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/bitrate.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/block.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/codebook.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/codebook.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/codec_internal.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/envelope.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/envelope.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/floor0.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/floor1.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/highlevel.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/info.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lookup.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lookup.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lookup_data.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lpc.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lpc.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lsp.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/lsp.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/mapping0.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/masking.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/mdct.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/mdct.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/misc.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/os.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/psy.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/psy.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/registry.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/registry.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/res0.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/scales.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/sharedbook.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/smallft.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/smallft.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/synthesis.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/vorbisenc.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/vorbisfile.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/window.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/libvorbis-1.3.2/lib/window.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/bitwise.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/codec.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/config_types.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/framing.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/Ogg Vorbis Licence.txt" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/ogg.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/os_types.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/vorbisenc.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/oggvorbis/vorbisfile.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_AiffAudioFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_AiffAudioFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_CoreAudioFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_CoreAudioFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_FlacAudioFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_FlacAudioFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_LAMEEncoderAudioFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_LAMEEncoderAudioFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_MP3AudioFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_MP3AudioFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_OggVorbisAudioFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_OggVorbisAudioFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_WavAudioFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_WavAudioFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_WindowsMediaAudioFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/codecs/juce_WindowsMediaAudioFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatManager.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatManager.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatReader.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatReader.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatReaderSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatReaderSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatWriter.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioFormatWriter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioSubsectionReader.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_AudioSubsectionReader.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_BufferingAudioFormatReader.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_BufferingAudioFormatReader.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/format/juce_MemoryMappedAudioFormatReader.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/sampler/juce_Sampler.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/sampler/juce_Sampler.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/juce_audio_formats.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/juce_audio_formats.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_formats/juce_audio_formats.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AAX/juce_AAX_Modifier_Injector.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AAX/juce_AAX_Wrapper.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUBase.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUBaseHelper.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUBuffer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUBuffer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUCarbonViewBase.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUCarbonViewBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUCarbonViewControl.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUCarbonViewControl.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUCarbonViewDispatch.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUDispatch.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUDispatch.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUInputElement.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUInputElement.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUInputFormatConverter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUMIDIBase.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUMIDIBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUOutputBase.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUOutputBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUOutputElement.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUOutputElement.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUPlugInDispatch.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUPlugInDispatch.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUScopeElement.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUScopeElement.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUSilentTimeout.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUTimestampGenerator.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/AUViewLocalizedStringKeys.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAtomic.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAtomicStack.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAudioChannelLayout.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAudioChannelLayout.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAUParameter.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAUParameter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAAutoDisposer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CADebugMacros.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CADebugPrintf.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAException.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAHostTimeBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAMath.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAMutex.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAMutex.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CarbonEventHandler.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CarbonEventHandler.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAReferenceCounted.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAStreamBasicDescription.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAStreamBasicDescription.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAThreadSafeList.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAVectorUnit.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAVectorUnit.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAVectorUnitTypes.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/CAXException.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/ComponentBase.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/ComponentBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/MusicDeviceBase.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/CoreAudioUtilityClasses/MusicDeviceBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/juce_AU_Wrapper.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/AU/juce_AUv3_Wrapper.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_DigiCode_Header.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_DigiCode1.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_DigiCode2.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_DigiCode3.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_MacUtilities.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_WinUtilities.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/RTAS/juce_RTAS_Wrapper.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/Standalone/juce_StandaloneFilterApp.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/Standalone/juce_StandaloneFilterWindow.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/Unity/juce_Unity_Wrapper.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/Unity/juce_UnityPluginInterface.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_CarbonVisibility.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_CheckSettingMacros.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_FakeMouseMoveGenerator.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_IncludeModuleHeaders.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_IncludeSystemHeaders.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_PluginHostType.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_PluginUtilities.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/utility/juce_WindowsHooks.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/VST/juce_VST_Wrapper.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/VST/juce_VST_Wrapper.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/VST/juce_VSTCallbackHandler.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/VST3/juce_VST3_Wrapper.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_AAX.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_AAX.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_AU_1.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_AU_2.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_AUv3.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_1.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_2.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_3.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_4.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_utils.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_RTAS_utils.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_Standalone.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_Unity.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_utils.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_VST_utils.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_VST2.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client_VST3.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_plugin_client/juce_audio_plugin_client.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format/juce_AudioPluginFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format/juce_AudioPluginFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format/juce_AudioPluginFormatManager.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format/juce_AudioPluginFormatManager.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/baseiids.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/classfactoryhelpers.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fbuffer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fbuffer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fdebug.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fdebug.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fobject.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fobject.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fstreamer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fstreamer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fstring.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/fstring.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/updatehandler.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/source/updatehandler.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/thread/include/flock.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/thread/source/flock.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/base/LICENSE.txt" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/conststringtable.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/conststringtable.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/coreiids.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/falignpop.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/falignpush.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/fplatform.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/fstrdefs.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ftypes.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/funknown.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/funknown.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/futils.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/fvariant.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ibstream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/icloneable.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ipersistent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ipluginbase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/istringresult.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/iupdatehandler.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/smartpointer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ustring.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/base/ustring.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/gui/iplugview.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/gui/iplugviewcontentscalesupport.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstattributes.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstaudioprocessor.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstautomationstate.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstchannelcontextinfo.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstcomponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstcontextmenu.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivsteditcontroller.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstevents.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivsthostapplication.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstinterappaudio.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstmessage.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstmidicontrollers.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstnoteexpression.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstparameterchanges.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstplugview.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstprefetchablesupport.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstprocesscontext.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstrepresentation.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/ivstunits.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/vstpshpack4.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/vstspeaker.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/vst/vsttypes.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/pluginterfaces/LICENSE.txt" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/common/memorystream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/common/memorystream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/common/pluginview.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/common/pluginview.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/hosting/hostclasses.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/hosting/hostclasses.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstbus.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstbus.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstcomponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstcomponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstcomponentbase.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstcomponentbase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vsteditcontroller.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vsteditcontroller.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstinitiids.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstparameters.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstparameters.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstpresetfile.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/source/vst/vstpresetfile.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/public.sdk/LICENSE.txt" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/LICENSE.txt" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/VST3_SDK/README.md" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_AU_Shared.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_AudioUnitPluginFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_AudioUnitPluginFormat.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_LADSPAPluginFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_LADSPAPluginFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_LegacyAudioParameter.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VST3Common.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VST3Headers.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VST3PluginFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VST3PluginFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VSTCommon.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VSTMidiEventList.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VSTPluginFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/format_types/juce_VSTPluginFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioPluginInstance.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioPluginInstance.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessor.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessor.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorEditor.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorEditor.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorGraph.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorGraph.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorListener.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorParameter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorParameterGroup.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_AudioProcessorParameterGroup.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_GenericAudioProcessorEditor.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_GenericAudioProcessorEditor.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_PluginDescription.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/processors/juce_PluginDescription.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_KnownPluginList.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_KnownPluginList.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_PluginDirectoryScanner.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_PluginDirectoryScanner.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_PluginListComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/scanning/juce_PluginListComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioParameterBool.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioParameterChoice.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioParameterFloat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioParameterInt.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioProcessorParameters.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioProcessorParameterWithID.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioProcessorValueTreeState.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_AudioProcessorValueTreeState.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/utilities/juce_RangedAudioParameter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/juce_audio_processors.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/juce_audio_processors.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_processors/juce_audio_processors.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/audio_cd/juce_AudioCDBurner.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/audio_cd/juce_AudioCDReader.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/audio_cd/juce_AudioCDReader.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioAppComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioAppComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioDeviceSelectorComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioDeviceSelectorComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioThumbnail.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioThumbnail.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioThumbnailBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioThumbnailCache.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioThumbnailCache.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioVisualiserComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_AudioVisualiserComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_BluetoothMidiDevicePairingDialogue.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_MidiKeyboardComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/gui/juce_MidiKeyboardComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_android_BluetoothMidiDevicePairingDialogue.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_ios_BluetoothMidiDevicePairingDialogue.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_linux_AudioCDReader.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_linux_BluetoothMidiDevicePairingDialogue.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_mac_AudioCDBurner.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_mac_AudioCDReader.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_mac_BluetoothMidiDevicePairingDialogue.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_win_BluetoothMidiDevicePairingDialogue.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_win32_AudioCDBurner.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/native/juce_win32_AudioCDReader.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/players/juce_AudioProcessorPlayer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/players/juce_AudioProcessorPlayer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/players/juce_SoundPlayer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/players/juce_SoundPlayer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/juce_audio_utils.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/juce_audio_utils.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_audio_utils/juce_audio_utils.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_AbstractFifo.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_AbstractFifo.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_Array.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ArrayAllocationBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ArrayBase.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ArrayBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_DynamicObject.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_DynamicObject.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ElementComparator.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_HashMap.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_HashMap_test.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_LinkedListPointer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ListenerList.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_NamedValueSet.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_NamedValueSet.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_OwnedArray.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_OwnedArray.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_PropertySet.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_PropertySet.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ReferenceCountedArray.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ReferenceCountedArray.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_ScopedValueSetter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_SortedSet.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_SparseSet.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_SparseSet.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_Variant.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/containers/juce_Variant.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_DirectoryIterator.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_DirectoryIterator.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_File.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_File.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileFilter.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileFilter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileInputStream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileInputStream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileOutputStream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileOutputStream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileSearchPath.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_FileSearchPath.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_MemoryMappedFile.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_TemporaryFile.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_TemporaryFile.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_WildcardFileFilter.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/files/juce_WildcardFileFilter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/javascript/juce_Javascript.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/javascript/juce_Javascript.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/javascript/juce_JSON.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/javascript/juce_JSON.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/logging/juce_FileLogger.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/logging/juce_FileLogger.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/logging/juce_Logger.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/logging/juce_Logger.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_BigInteger.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_BigInteger.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_Expression.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_Expression.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_MathsFunctions.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_NormalisableRange.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_Random.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_Random.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_Range.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/maths/juce_StatisticsAccumulator.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_Atomic.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_ByteOrder.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_ContainerDeletePolicy.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_HeapBlock.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_HeavyweightLeakedObjectDetector.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_LeakedObjectDetector.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_Memory.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_MemoryBlock.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_MemoryBlock.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_OptionalScopedPointer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_ReferenceCountedObject.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_ScopedPointer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_SharedResourcePointer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_Singleton.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/memory/juce_WeakReference.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_ConsoleApplication.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_ConsoleApplication.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_Result.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_Result.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_RuntimePermissions.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_RuntimePermissions.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_StdFunctionCompat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_StdFunctionCompat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_Uuid.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_Uuid.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/misc/juce_WindowsRegistry.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_Files.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_JNIHelpers.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_Misc.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_Network.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_RuntimePermissions.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_SystemStats.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_android_Threads.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_BasicNativeHeaders.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_curl_Network.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_linux_CommonFile.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_linux_Files.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_linux_Network.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_linux_SystemStats.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_linux_Threads.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_ClangBugWorkaround.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_Files.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_Network.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_Strings.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_SystemStats.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_mac_Threads.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_osx_ObjCHelpers.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_posix_IPAddress.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_posix_NamedPipe.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_posix_SharedCode.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_ComSmartPtr.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_Files.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_Network.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_Registry.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_SystemStats.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/native/juce_win32_Threads.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_IPAddress.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_IPAddress.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_MACAddress.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_MACAddress.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_NamedPipe.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_NamedPipe.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_Socket.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_Socket.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_URL.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_URL.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_WebInputStream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/network/juce_WebInputStream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_BufferedInputStream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_BufferedInputStream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_FileInputSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_FileInputSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_InputSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_InputStream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_InputStream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_MemoryInputStream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_MemoryInputStream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_MemoryOutputStream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_MemoryOutputStream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_OutputStream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_OutputStream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_SubregionStream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_SubregionStream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_URLInputSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/streams/juce_URLInputSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_CompilerSupport.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_PlatformDefs.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_StandardHeader.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_SystemStats.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_SystemStats.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/system/juce_TargetPlatform.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_Base64.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_Base64.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharacterFunctions.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharacterFunctions.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharPointer_ASCII.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharPointer_UTF8.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharPointer_UTF16.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_CharPointer_UTF32.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_Identifier.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_Identifier.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_LocalisedStrings.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_LocalisedStrings.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_NewLine.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_String.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_String.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringArray.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringArray.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringPairArray.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringPairArray.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringPool.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringPool.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_StringRef.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_TextDiff.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/text/juce_TextDiff.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ChildProcess.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ChildProcess.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_CriticalSection.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_DynamicLibrary.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_HighResolutionTimer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_HighResolutionTimer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_InterProcessLock.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_Process.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ReadWriteLock.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ReadWriteLock.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ScopedLock.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ScopedReadLock.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ScopedWriteLock.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_SpinLock.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_Thread.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_Thread.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ThreadLocalValue.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ThreadPool.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_ThreadPool.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_TimeSliceThread.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_TimeSliceThread.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/threads/juce_WaitableEvent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_PerformanceCounter.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_PerformanceCounter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_RelativeTime.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_RelativeTime.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_Time.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/time/juce_Time.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/unit_tests/juce_UnitTest.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/unit_tests/juce_UnitTest.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/xml/juce_XmlDocument.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/xml/juce_XmlDocument.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/xml/juce_XmlElement.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/xml/juce_XmlElement.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/adler32.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/compress.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/crc32.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/crc32.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/deflate.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/deflate.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/infback.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inffast.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inffast.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inffixed.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inflate.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inflate.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inftrees.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/inftrees.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/trees.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/trees.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/uncompr.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/zconf.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/zconf.in.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/zlib.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/zutil.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/zlib/zutil.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_GZIPCompressorOutputStream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_GZIPCompressorOutputStream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_GZIPDecompressorInputStream.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_GZIPDecompressorInputStream.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_ZipFile.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/zip/juce_ZipFile.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/juce_core.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/juce_core.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_core/juce_core.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/app_properties/juce_ApplicationProperties.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/app_properties/juce_ApplicationProperties.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/app_properties/juce_PropertiesFile.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/app_properties/juce_PropertiesFile.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/undomanager/juce_UndoableAction.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/undomanager/juce_UndoManager.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/undomanager/juce_UndoManager.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_CachedValue.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_CachedValue.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_Value.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_Value.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_ValueTree.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_ValueTree.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_ValueTreeSynchroniser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_ValueTreeSynchroniser.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/values/juce_ValueWithDefault.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/juce_data_structures.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/juce_data_structures.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_data_structures/juce_data_structures.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ActionBroadcaster.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ActionBroadcaster.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ActionListener.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_AsyncUpdater.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_AsyncUpdater.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ChangeBroadcaster.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ChangeBroadcaster.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/broadcasters/juce_ChangeListener.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_ConnectedChildProcess.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_ConnectedChildProcess.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_InterprocessConnection.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_InterprocessConnection.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_InterprocessConnectionServer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_InterprocessConnectionServer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_NetworkServiceDiscovery.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/interprocess/juce_NetworkServiceDiscovery.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_ApplicationBase.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_ApplicationBase.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_CallbackMessage.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_DeletedAtShutdown.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_DeletedAtShutdown.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_Initialisation.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_Message.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_MessageListener.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_MessageListener.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_MessageManager.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_MessageManager.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_MountedVolumeListChangeDetector.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/messages/juce_NotificationType.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_android_Messaging.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_ios_MessageManager.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_linux_EventLoop.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_linux_Messaging.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_mac_MessageManager.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_osx_MessageQueue.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_win32_HiddenMessageWindow.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_win32_Messaging.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_win32_WinRTWrapper.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/native/juce_win32_WinRTWrapper.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/timers/juce_MultiTimer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/timers/juce_MultiTimer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/timers/juce_Timer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/timers/juce_Timer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/juce_events.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/juce_events.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_events/juce_events.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_Colour.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_Colour.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_ColourGradient.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_ColourGradient.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_Colours.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_Colours.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_FillType.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_FillType.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/colour/juce_PixelFormats.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_GraphicsContext.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_GraphicsContext.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_LowLevelGraphicsContext.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_LowLevelGraphicsPostScriptRenderer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_LowLevelGraphicsPostScriptRenderer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_LowLevelGraphicsSoftwareRenderer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/contexts/juce_LowLevelGraphicsSoftwareRenderer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/effects/juce_DropShadowEffect.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/effects/juce_DropShadowEffect.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/effects/juce_GlowEffect.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/effects/juce_GlowEffect.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/effects/juce_ImageEffectFilter.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_AttributedString.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_AttributedString.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_CustomTypeface.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_CustomTypeface.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_Font.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_Font.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_GlyphArrangement.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_GlyphArrangement.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_TextLayout.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_TextLayout.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_Typeface.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/fonts/juce_Typeface.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_AffineTransform.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_AffineTransform.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_BorderSize.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_EdgeTable.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_EdgeTable.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Line.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Parallelogram.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Path.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Path.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_PathIterator.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_PathIterator.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_PathStrokeType.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_PathStrokeType.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Point.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_Rectangle.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/geometry/juce_RectangleList.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/cderror.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/changes to libjpeg for JUCE.txt" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcapimin.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcapistd.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jccoefct.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jccolor.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcdctmgr.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jchuff.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jchuff.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcinit.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcmainct.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcmarker.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcmaster.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcomapi.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jconfig.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcparam.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcphuff.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcprepct.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jcsample.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jctrans.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdapimin.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdapistd.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdatasrc.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdcoefct.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdcolor.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdct.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jddctmgr.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdhuff.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdhuff.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdinput.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdmainct.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdmarker.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdmaster.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdmerge.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdphuff.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdpostct.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdsample.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jdtrans.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jerror.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jerror.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jfdctflt.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jfdctfst.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jfdctint.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jidctflt.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jidctfst.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jidctint.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jidctred.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jinclude.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jmemmgr.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jmemnobs.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jmemsys.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jmorecfg.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jpegint.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jpeglib.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jquant1.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jquant2.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jutils.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/jversion.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/transupp.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/jpglib/transupp.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/libpng_readme.txt" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/png.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/png.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngconf.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngerror.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngget.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pnginfo.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngmem.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngpread.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngpriv.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngread.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngrio.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngrtran.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngrutil.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngset.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngstruct.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngtrans.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngwio.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngwrite.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngwtran.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/pnglib/pngwutil.c" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/juce_GIFLoader.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/juce_JPEGLoader.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/image_formats/juce_PNGLoader.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_Image.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_Image.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageCache.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageCache.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageConvolutionKernel.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageConvolutionKernel.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageFileFormat.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/images/juce_ImageFileFormat.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_android_Fonts.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_android_GraphicsContext.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_android_IconHelpers.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_freetype_Fonts.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_linux_Fonts.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_linux_IconHelpers.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_mac_CoreGraphicsContext.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_mac_CoreGraphicsContext.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_mac_CoreGraphicsHelpers.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_mac_Fonts.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_mac_IconHelpers.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_RenderingHelpers.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_Direct2DGraphicsContext.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_Direct2DGraphicsContext.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_DirectWriteTypeface.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_DirectWriteTypeLayout.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_Fonts.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/native/juce_win32_IconHelpers.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/placement/juce_Justification.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/placement/juce_RectanglePlacement.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/placement/juce_RectanglePlacement.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/juce_graphics.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/juce_graphics.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_graphics/juce_graphics.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/application/juce_Application.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/application/juce_Application.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ArrowButton.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ArrowButton.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_Button.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_Button.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_DrawableButton.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_DrawableButton.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_HyperlinkButton.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_HyperlinkButton.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ImageButton.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ImageButton.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ShapeButton.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ShapeButton.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_TextButton.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_TextButton.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ToggleButton.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ToggleButton.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ToolbarButton.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/buttons/juce_ToolbarButton.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandID.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandInfo.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandInfo.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandManager.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandManager.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandTarget.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_ApplicationCommandTarget.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_KeyPressMappingSet.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/commands/juce_KeyPressMappingSet.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_CachedComponentImage.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_Component.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_Component.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_ComponentListener.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_ComponentListener.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_ModalComponentManager.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/components/juce_ModalComponentManager.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/desktop/juce_Desktop.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/desktop/juce_Desktop.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/desktop/juce_Displays.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/desktop/juce_Displays.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_Drawable.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_Drawable.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableComposite.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableComposite.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableImage.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableImage.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawablePath.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawablePath.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableRectangle.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableRectangle.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableShape.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableShape.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableText.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_DrawableText.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/drawables/juce_SVGParser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_ContentSharer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_ContentSharer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_DirectoryContentsDisplayComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_DirectoryContentsDisplayComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_DirectoryContentsList.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_DirectoryContentsList.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileBrowserComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileBrowserComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileBrowserListener.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileChooser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileChooser.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileChooserDialogBox.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileChooserDialogBox.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileListComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileListComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FilenameComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FilenameComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FilePreviewComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileSearchPathListComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileSearchPathListComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileTreeComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_FileTreeComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_ImagePreviewComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/filebrowser/juce_ImagePreviewComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_CaretComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_CaretComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyboardFocusTraverser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyboardFocusTraverser.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyListener.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyListener.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyPress.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_KeyPress.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_ModifierKeys.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_ModifierKeys.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_SystemClipboard.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_TextEditorKeyMapper.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/keyboard/juce_TextInputTarget.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_AnimatedPosition.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_AnimatedPositionBehaviours.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentAnimator.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentAnimator.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentBoundsConstrainer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentBoundsConstrainer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentBuilder.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentBuilder.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentMovementWatcher.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ComponentMovementWatcher.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ConcertinaPanel.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ConcertinaPanel.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_FlexBox.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_FlexBox.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_FlexItem.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_Grid.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_Grid.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_GridItem.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_GridItem.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_GridUnitTests.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_GroupComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_GroupComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_MultiDocumentPanel.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_MultiDocumentPanel.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableBorderComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableBorderComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableCornerComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableCornerComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableEdgeComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ResizableEdgeComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ScrollBar.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_ScrollBar.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_SidePanel.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_SidePanel.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableLayoutManager.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableLayoutManager.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableLayoutResizerBar.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableLayoutResizerBar.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableObjectResizer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_StretchableObjectResizer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_TabbedButtonBar.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_TabbedButtonBar.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_TabbedComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_TabbedComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_Viewport.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/layout/juce_Viewport.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V1.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V1.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V2.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V2.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V3.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V3.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V4.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/lookandfeel/juce_LookAndFeel_V4.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_BurgerMenuComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_BurgerMenuComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_MenuBarComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_MenuBarComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_MenuBarModel.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_MenuBarModel.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_PopupMenu.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/menus/juce_PopupMenu.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_BubbleComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_BubbleComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_DropShadower.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_DropShadower.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_JUCESplashScreen.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/misc/juce_JUCESplashScreen.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_ComponentDragger.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_ComponentDragger.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_DragAndDropContainer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_DragAndDropContainer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_DragAndDropTarget.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_FileDragAndDropTarget.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_LassoComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseCursor.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseCursor.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseEvent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseEvent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseInactivityDetector.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseInactivityDetector.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseInputSource.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseInputSource.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseListener.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_MouseListener.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_SelectedItemSet.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_TextDragAndDropTarget.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/mouse/juce_TooltipClient.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_android_ContentSharer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_android_FileChooser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_android_Windowing.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_common_MimeTypes.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_ios_ContentSharer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_ios_FileChooser.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_ios_UIViewComponentPeer.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_ios_Windowing.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_linux_FileChooser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_linux_X11.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_linux_X11.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_linux_X11_Clipboard.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_linux_X11_Windowing.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_mac_FileChooser.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_mac_MainMenu.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_mac_MouseCursor.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_mac_NSViewComponentPeer.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_mac_Windowing.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_MultiTouchMapper.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_win32_DragAndDrop.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_win32_FileChooser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/native/juce_win32_Windowing.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_MarkerList.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_MarkerList.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeCoordinate.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeCoordinate.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeCoordinatePositioner.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeCoordinatePositioner.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeParallelogram.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeParallelogram.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativePoint.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativePoint.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativePointPath.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativePointPath.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeRectangle.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/positioning/juce_RelativeRectangle.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_BooleanPropertyComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_BooleanPropertyComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_ButtonPropertyComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_ButtonPropertyComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_ChoicePropertyComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_ChoicePropertyComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_MultiChoicePropertyComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_MultiChoicePropertyComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_PropertyComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_PropertyComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_PropertyPanel.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_PropertyPanel.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_SliderPropertyComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_SliderPropertyComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_TextPropertyComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/properties/juce_TextPropertyComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ComboBox.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ComboBox.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ImageComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ImageComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Label.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Label.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ListBox.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ListBox.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ProgressBar.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ProgressBar.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Slider.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Slider.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TableHeaderComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TableHeaderComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TableListBox.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TableListBox.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TextEditor.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TextEditor.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Toolbar.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_Toolbar.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ToolbarItemComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ToolbarItemComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ToolbarItemFactory.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ToolbarItemPalette.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_ToolbarItemPalette.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TreeView.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/widgets/juce_TreeView.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_AlertWindow.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_AlertWindow.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_CallOutBox.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_CallOutBox.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ComponentPeer.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ComponentPeer.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_DialogWindow.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_DialogWindow.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_DocumentWindow.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_DocumentWindow.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_NativeMessageBox.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ResizableWindow.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ResizableWindow.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ThreadWithProgressWindow.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_ThreadWithProgressWindow.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_TooltipWindow.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_TooltipWindow.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_TopLevelWindow.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/windows/juce_TopLevelWindow.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/juce_gui_basics.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/juce_gui_basics.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_basics/juce_gui_basics.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CodeDocument.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CodeDocument.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CodeEditorComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CodeEditorComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CodeTokeniser.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CPlusPlusCodeTokeniser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CPlusPlusCodeTokeniser.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_CPlusPlusCodeTokeniserFunctions.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_LuaCodeTokeniser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_LuaCodeTokeniser.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_XMLCodeTokeniser.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/code_editor/juce_XMLCodeTokeniser.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/documents/juce_FileBasedDocument.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/documents/juce_FileBasedDocument.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_ActiveXControlComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_AndroidViewComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_NSViewComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_ScopedDPIAwarenessDisabler.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_UIViewComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/embedding/juce_XEmbedComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_AnimatedAppComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_AnimatedAppComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_AppleRemote.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_BubbleMessageComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_BubbleMessageComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_ColourSelector.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_ColourSelector.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_KeyMappingEditorComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_KeyMappingEditorComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_LiveConstantEditor.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_LiveConstantEditor.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_PreferencesPanel.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_PreferencesPanel.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_PushNotifications.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_PushNotifications.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_RecentlyOpenedFilesList.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_RecentlyOpenedFilesList.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_SplashScreen.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_SplashScreen.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_SystemTrayIconComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_SystemTrayIconComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/misc/juce_WebBrowserComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_android_PushNotifications.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_android_WebBrowserComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_AndroidViewComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_ios_PushNotifications.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_ios_UIViewComponent.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_linux_X11_SystemTrayIcon.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_linux_X11_WebBrowserComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_linux_XEmbedComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_AppleRemote.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_CarbonViewWrapperComponent.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_NSViewComponent.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_PushNotifications.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_SystemTrayIcon.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_mac_WebBrowserComponent.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_win32_ActiveXComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_win32_SystemTrayIcon.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/native/juce_win32_WebBrowserComponent.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/juce_gui_extra.cpp" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/juce_gui_extra.mm" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../../../../newJUCE/JUCE/modules/juce_gui_extra/juce_gui_extra.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../JuceLibraryCode/AppConfig.h" PROPERTIES HEADER_FILE_ONLY TRUE)
set_source_files_properties("../../../JuceLibraryCode/JuceHeader.h" PROPERTIES HEADER_FILE_ONLY TRUE)

find_library(log "log")
find_library(android "android")
find_library(glesv2 "GLESv2")
find_library(egl "EGL")

target_link_libraries( ${BINARY_NAME}

    ${log}
    ${android}
    ${glesv2}
    ${egl}
    "cpufeatures"
)
