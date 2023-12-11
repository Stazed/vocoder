// DPF includes
#include "DistrhoPlugin.hpp"
#include "DSP.hpp"

START_NAMESPACE_DISTRHO

class Vocoder : public Plugin {
    public:
        Vocoder() : Plugin(kParameterCount, 0, 0), bypass(0.0)
        {
            v_process = new VocProc(Plugin::getSampleRate());
        }
        
        ~Vocoder()
        {
            if(v_process)
                delete v_process;
        }
    protected:

        const char *getLabel() const override { return "Vocoder"; }
        const char *getDescription() const override {
            return "Simple vocoder plugin.";
        }
        const char *getMaker() const override { return "Stazed"; }
        const char *getLicense() const override { return "GPL2"; }
        uint32_t getVersion() const override { return d_version(1,0,0); }
        int64_t getUniqueId() const override { 
            return d_cconst('V','O','C','T'); 
        }
    
    void initParameter (uint32_t index, Parameter& parameter) override {
        switch (index) {
            case kBypass:
                parameter.hints = kParameterIsAutomatable|kParameterIsBoolean|kParameterIsInteger;
                parameter.name = "Bypass";
                parameter.symbol = "bypass";
                parameter.ranges.def = 0.0f;
                parameter.ranges.min = 0.0f;
                parameter.ranges.max = 1.0f;
                break;
            default:
                break;
        }
    }

    float getParameterValue(uint32_t index) const override {
        switch (index) {
        case kBypass:
            return bypass;
        default:
            return 0.0;
        }
    }

    void setParameterValue(uint32_t index, float value) override
    {
        switch (index)
        {
        case kBypass:
            bypass = value;
            v_process->set_bypass(bypass);
            break;
        default:
            break;
        }
    }

    void run(const float **inputs, float **outputs, uint32_t frames) override
    {
        v_process->run(inputs, outputs, frames);
    }

    private:
        float bypass;
        VocProc *v_process;

        DISTRHO_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(Vocoder);
};


Plugin *createPlugin()
{
    return new Vocoder();
}

END_NAMESPACE_DISTRHO
