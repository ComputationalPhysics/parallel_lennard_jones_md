#pragma once
class System;
class StatisticsSampler;

class Thermostat
{
public:
    double relaxation_time;

    Thermostat(double relaxation_time_);
    void apply(StatisticsSampler *sampler, System *system, const double &temperature, bool frozen_atoms_only);
};
