// DPF includes
#include "DistrhoPlugin.hpp"
START_NAMESPACE_DISTRHO

class voc : public Plugin {
    public:
        voc() : Plugin(kParameterCount, 0, 0), bypass(1.0) {}
    protected:

        const char *getLabel() const override { return "voc"; }
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
                parameter.name = "Bypass";
                parameter.symbol = "bypass";
                parameter.ranges.def = 1.0f;
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
            return 1.0;
        }
    }

    void setParameterValue(uint32_t index, float value) override {
        switch (index) {
        case kBypass:
            bypass = value;
            break;
        default:
            break;
        }
    }

    void run(const float **inputs, float **outputs, uint32_t frames) override {
        const float *const in = inputs[0];
        float *const out = outputs[0];

        for (uint32_t i = 0; i < frames; i++) {
            out[i] = in[i] * bypass;
        }
    }

    private:
        float bypass;

        DISTRHO_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(voc);
};


Plugin *createPlugin()
{
    return new voc();
}

END_NAMESPACE_DISTRHO
