#ifdef MEMORY_POOL
#endif /* MEMORY_POOL */ 
#ifdef DEFINE_STORAGE
#endif /* DEFINE_STORAGE */ 
#ifdef NEW_PROBLEM_TYPES
#endif /* NEW_PROBLEM_TYPES */ 

    if (strncmp(parameter_name, "NumberOfGhostZones", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfGhostZones);
        return rv;
    }

    if (strncmp(parameter_name, "LoadBalancing", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) LoadBalancing);
        return rv;
    }

    if (strncmp(parameter_name, "LoadBalancingCycleSkip", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) LoadBalancingCycleSkip);
        return rv;
    }

    if (strncmp(parameter_name, "ResetLoadBalancing", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ResetLoadBalancing);
        return rv;
    }

    if (strncmp(parameter_name, "CoresPerNode", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) CoresPerNode);
        return rv;
    }

    if (strncmp(parameter_name, "PreviousMaxTask", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) PreviousMaxTask);
        return rv;
    }

    if (strncmp(parameter_name, "LoadBalancingMinLevel", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) LoadBalancingMinLevel);
        return rv;
    }

    if (strncmp(parameter_name, "LoadBalancingMaxLevel", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) LoadBalancingMaxLevel);
        return rv;
    }

    if (strncmp(parameter_name, "FileDirectedOutput", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) FileDirectedOutput);
        return rv;
    }

    if (strncmp(parameter_name, "HierarchyFileInputFormat", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) HierarchyFileInputFormat);
        return rv;
    }

    if (strncmp(parameter_name, "HierarchyFileOutputFormat", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) HierarchyFileOutputFormat);
        return rv;
    }

    if (strncmp(parameter_name, "TotalNumberOfGrids", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) TotalNumberOfGrids);
        return rv;
    }

    if (strncmp(parameter_name, "debug", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) debug);
        return rv;
    }

    if (strncmp(parameter_name, "debug1", 7) == 0) {
        PyObject *rv = PyLong_FromLong((long) debug1);
        return rv;
    }

    if (strncmp(parameter_name, "debug2", 7) == 0) {
        PyObject *rv = PyLong_FromLong((long) debug2);
        return rv;
    }

    if (strncmp(parameter_name, "extract", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) extract);
        return rv;
    }

    if (strncmp(parameter_name, "CheckpointRestart", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) CheckpointRestart);
        return rv;
    }

    if (strncmp(parameter_name, "WriteGhostZones", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) WriteGhostZones);
        return rv;
    }

    if (strncmp(parameter_name, "ReadGhostZones", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) ReadGhostZones);
        return rv;
    }

    if (strncmp(parameter_name, "ProblemType", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) ProblemType);
        return rv;
    }
#ifdef NEW_PROBLEM_TYPES
#endif /* NEW_PROBLEM_TYPES */ 

    if (strncmp(parameter_name, "huge_number", 12) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) huge_number);
        return rv;
    }

    if (strncmp(parameter_name, "tiny_number", 12) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) tiny_number);
        return rv;
    }

    if (strncmp(parameter_name, "Gamma", 6) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) Gamma);
        return rv;
    }

    if (strncmp(parameter_name, "PressureFree", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) PressureFree);
        return rv;
    }

    if (strncmp(parameter_name, "QuantumPressure", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) QuantumPressure);
        return rv;
    }

    if (strncmp(parameter_name, "FDMMass", 8) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) FDMMass);
        return rv;
    }

    if (strncmp(parameter_name, "RefineBy", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) RefineBy);
        return rv;
    }

    if (strncmp(parameter_name, "MaximumRefinementLevel", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) MaximumRefinementLevel);
        return rv;
    }

    if (strncmp(parameter_name, "MaximumGravityRefinementLevel", 30) == 0) {
        PyObject *rv = PyLong_FromLong((long) MaximumGravityRefinementLevel);
        return rv;
    }

    if (strncmp(parameter_name, "MaximumParticleRefinementLevel", 31) == 0) {
        PyObject *rv = PyLong_FromLong((long) MaximumParticleRefinementLevel);
        return rv;
    }

    if (strncmp(parameter_name, "FastSiblingLocatorEntireDomain", 31) == 0) {
        PyObject *rv = PyLong_FromLong((long) FastSiblingLocatorEntireDomain);
        return rv;
    }

    if (strncmp(parameter_name, "MustRefineRegionMinRefinementLevel", 35) == 0) {
        PyObject *rv = PyLong_FromLong((long) MustRefineRegionMinRefinementLevel);
        return rv;
    }

    if (strncmp(parameter_name, "MetallicityRefinementMinLevel", 30) == 0) {
        PyObject *rv = PyLong_FromLong((long) MetallicityRefinementMinLevel);
        return rv;
    }

    if (strncmp(parameter_name, "MetallicityRefinementMinMetallicity", 36) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MetallicityRefinementMinMetallicity);
        return rv;
    }

    if (strncmp(parameter_name, "MetallicityRefinementMinDensity", 32) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MetallicityRefinementMinDensity);
        return rv;
    }

    if (strncmp(parameter_name, "TimestepSafetyVelocity", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) TimestepSafetyVelocity);
        return rv;
    }

    if (strncmp(parameter_name, "FluxCorrection", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) FluxCorrection);
        return rv;
    }

    if (strncmp(parameter_name, "UseCoolingTimestep", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseCoolingTimestep);
        return rv;
    }

    if (strncmp(parameter_name, "CoolingTimestepSafetyFactor", 28) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CoolingTimestepSafetyFactor);
        return rv;
    }

    if (strncmp(parameter_name, "ConservativeInterpolation", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) ConservativeInterpolation);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumEfficiency", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MinimumEfficiency);
        return rv;
    }

    if (strncmp(parameter_name, "SubgridSizeAutoAdjust", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) SubgridSizeAutoAdjust);
        return rv;
    }

    if (strncmp(parameter_name, "OptimalSubgridsPerProcessor", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) OptimalSubgridsPerProcessor);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumSubgridEdge", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) MinimumSubgridEdge);
        return rv;
    }

    if (strncmp(parameter_name, "MaximumSubgridSize", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) MaximumSubgridSize);
        return rv;
    }

    if (strncmp(parameter_name, "CriticalGridRatio", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CriticalGridRatio);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfBufferZones", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfBufferZones);
        return rv;
    }

    if (strncmp(parameter_name, "DimUnits", 9) == 0) {
        PyObject *rv = PyString_FromString((char *) DimUnits);
        return rv;
    }

    if (strncmp(parameter_name, "DimLabels", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) DimLabels);
        return rv;
    }

    if (strncmp(parameter_name, "DataLabel", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) DataLabel);
        return rv;
    }

    if (strncmp(parameter_name, "DataUnits", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) DataUnits);
        return rv;
    }

    if (strncmp(parameter_name, "RefineRegionAutoAdjust", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) RefineRegionAutoAdjust);
        return rv;
    }

    if (strncmp(parameter_name, "MultiRefineRegion", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) MultiRefineRegion);
        return rv;
    }

    if (strncmp(parameter_name, "MultiRefineRegionMaximumOuterLevel", 35) == 0) {
        PyObject *rv = PyLong_FromLong((long) MultiRefineRegionMaximumOuterLevel);
        return rv;
    }

    if (strncmp(parameter_name, "MultiRefineRegionMinimumOuterLevel", 35) == 0) {
        PyObject *rv = PyLong_FromLong((long) MultiRefineRegionMinimumOuterLevel);
        return rv;
    }

    if (strncmp(parameter_name, "UniformGravity", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) UniformGravity);
        return rv;
    }

    if (strncmp(parameter_name, "UniformGravityDirection", 24) == 0) {
        PyObject *rv = PyLong_FromLong((long) UniformGravityDirection);
        return rv;
    }

    if (strncmp(parameter_name, "UniformGravityConstant", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) UniformGravityConstant);
        return rv;
    }

    if (strncmp(parameter_name, "PointSourceGravity", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) PointSourceGravity);
        return rv;
    }

    if (strncmp(parameter_name, "PointSourceGravityConstant", 27) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) PointSourceGravityConstant);
        return rv;
    }

    if (strncmp(parameter_name, "PointSourceGravityCoreRadius", 29) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) PointSourceGravityCoreRadius);
        return rv;
    }

    if (strncmp(parameter_name, "DiskGravity", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) DiskGravity);
        return rv;
    }

    if (strncmp(parameter_name, "DiskGravityStellarDiskMass", 27) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DiskGravityStellarDiskMass);
        return rv;
    }

    if (strncmp(parameter_name, "DiskGravityStellarDiskScaleHeightR", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DiskGravityStellarDiskScaleHeightR);
        return rv;
    }

    if (strncmp(parameter_name, "DiskGravityStellarDiskScaleHeightz", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DiskGravityStellarDiskScaleHeightz);
        return rv;
    }

    if (strncmp(parameter_name, "DiskGravityStellarBulgeMass", 28) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DiskGravityStellarBulgeMass);
        return rv;
    }

    if (strncmp(parameter_name, "DiskGravityStellarBulgeR", 25) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DiskGravityStellarBulgeR);
        return rv;
    }

    if (strncmp(parameter_name, "DiskGravityDarkMatterR", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DiskGravityDarkMatterR);
        return rv;
    }

    if (strncmp(parameter_name, "DiskGravityDarkMatterDensity", 29) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DiskGravityDarkMatterDensity);
        return rv;
    }

    if (strncmp(parameter_name, "SelfGravity", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) SelfGravity);
        return rv;
    }

    if (strncmp(parameter_name, "SelfGravityGasOff", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) SelfGravityGasOff);
        return rv;
    }

    if (strncmp(parameter_name, "AccretionKernal", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) AccretionKernal);
        return rv;
    }

    if (strncmp(parameter_name, "CopyGravPotential", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) CopyGravPotential);
        return rv;
    }

    if (strncmp(parameter_name, "PotentialIterations", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) PotentialIterations);
        return rv;
    }

    if (strncmp(parameter_name, "BaryonSelfGravityApproximation", 31) == 0) {
        PyObject *rv = PyLong_FromLong((long) BaryonSelfGravityApproximation);
        return rv;
    }

    if (strncmp(parameter_name, "GravitationalConstant", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GravitationalConstant);
        return rv;
    }

    if (strncmp(parameter_name, "ComputePotential", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) ComputePotential);
        return rv;
    }

    if (strncmp(parameter_name, "WritePotential", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) WritePotential);
        return rv;
    }

    if (strncmp(parameter_name, "ParticleSubgridDepositMode", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) ParticleSubgridDepositMode);
        return rv;
    }

    if (strncmp(parameter_name, "DualEnergyFormalism", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) DualEnergyFormalism);
        return rv;
    }

    if (strncmp(parameter_name, "DualEnergyFormalismEta1", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DualEnergyFormalismEta1);
        return rv;
    }

    if (strncmp(parameter_name, "DualEnergyFormalismEta2", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DualEnergyFormalismEta2);
        return rv;
    }

    if (strncmp(parameter_name, "ParticleCourantSafetyNumber", 28) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ParticleCourantSafetyNumber);
        return rv;
    }

    if (strncmp(parameter_name, "RootGridCourantSafetyNumber", 28) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RootGridCourantSafetyNumber);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeCooling", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeCooling);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeCoolingModel", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeCoolingModel);
        return rv;
    }

    if (strncmp(parameter_name, "use_grackle", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) use_grackle);
        return rv;
    }

    if (strncmp(parameter_name, "GadgetEquilibriumCooling", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) GadgetEquilibriumCooling);
        return rv;
    }

    if (strncmp(parameter_name, "RandomForcing", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) RandomForcing);
        return rv;
    }

    if (strncmp(parameter_name, "RandomForcingEdot", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RandomForcingEdot);
        return rv;
    }

    if (strncmp(parameter_name, "RandomForcingMachNumber", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RandomForcingMachNumber);
        return rv;
    }

    if (strncmp(parameter_name, "DrivenFlowSeed", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) DrivenFlowSeed);
        return rv;
    }

    if (strncmp(parameter_name, "DrivenFlowWeight", 17) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DrivenFlowWeight);
        return rv;
    }

    if (strncmp(parameter_name, "UseSGSModel", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseSGSModel);
        return rv;
    }

    if (strncmp(parameter_name, "SGSFilterStencil", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) SGSFilterStencil);
        return rv;
    }

    if (strncmp(parameter_name, "SGSNeedJacobians", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) SGSNeedJacobians);
        return rv;
    }

    if (strncmp(parameter_name, "SGSNeedMixedFilteredQuantities", 31) == 0) {
        PyObject *rv = PyLong_FromLong((long) SGSNeedMixedFilteredQuantities);
        return rv;
    }

    if (strncmp(parameter_name, "SGSFilterWidth", 15) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SGSFilterWidth);
        return rv;
    }

    if (strncmp(parameter_name, "SGScoeffERS2M2Star", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SGScoeffERS2M2Star);
        return rv;
    }

    if (strncmp(parameter_name, "SGScoeffEVStarEnS2Star", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SGScoeffEVStarEnS2Star);
        return rv;
    }

    if (strncmp(parameter_name, "SGScoeffEnS2StarTrace", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SGScoeffEnS2StarTrace);
        return rv;
    }

    if (strncmp(parameter_name, "SGScoeffNLemfCompr", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SGScoeffNLemfCompr);
        return rv;
    }

    if (strncmp(parameter_name, "SGScoeffNLu", 12) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SGScoeffNLu);
        return rv;
    }

    if (strncmp(parameter_name, "SGScoeffNLuNormedEnS2Star", 26) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SGScoeffNLuNormedEnS2Star);
        return rv;
    }

    if (strncmp(parameter_name, "SGScoeffNLb", 12) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SGScoeffNLb);
        return rv;
    }

    if (strncmp(parameter_name, "SGScoeffSSu", 12) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SGScoeffSSu);
        return rv;
    }

    if (strncmp(parameter_name, "SGScoeffSSb", 12) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SGScoeffSSb);
        return rv;
    }

    if (strncmp(parameter_name, "SGScoeffSSemf", 14) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SGScoeffSSemf);
        return rv;
    }

    if (strncmp(parameter_name, "MultiSpecies", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) MultiSpecies);
        return rv;
    }

    if (strncmp(parameter_name, "NoMultiSpeciesButColors", 24) == 0) {
        PyObject *rv = PyLong_FromLong((long) NoMultiSpeciesButColors);
        return rv;
    }

    if (strncmp(parameter_name, "ThreeBodyRate", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) ThreeBodyRate);
        return rv;
    }

    if (strncmp(parameter_name, "H2FormationOnDust", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) H2FormationOnDust);
        return rv;
    }

    if (strncmp(parameter_name, "MixSpeciesAndColors", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) MixSpeciesAndColors);
        return rv;
    }

    if (strncmp(parameter_name, "GloverChemistryModel", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) GloverChemistryModel);
        return rv;
    }

    if (strncmp(parameter_name, "MultiMetals", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) MultiMetals);
        return rv;
    }

    if (strncmp(parameter_name, "CRModel", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) CRModel);
        return rv;
    }

    if (strncmp(parameter_name, "CRDiffusion", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) CRDiffusion);
        return rv;
    }

    if (strncmp(parameter_name, "CRFeedback", 11) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CRFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "CRkappa", 8) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CRkappa);
        return rv;
    }

    if (strncmp(parameter_name, "CRCourantSafetyNumber", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CRCourantSafetyNumber);
        return rv;
    }

    if (strncmp(parameter_name, "CRdensFloor", 12) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CRdensFloor);
        return rv;
    }

    if (strncmp(parameter_name, "CRmaxSoundSpeed", 16) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CRmaxSoundSpeed);
        return rv;
    }

    if (strncmp(parameter_name, "CRgamma", 8) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CRgamma);
        return rv;
    }

    if (strncmp(parameter_name, "CosmologySimulationUniformCR", 29) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CosmologySimulationUniformCR);
        return rv;
    }

    if (strncmp(parameter_name, "ShockMethod", 12) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShockMethod);
        return rv;
    }

    if (strncmp(parameter_name, "ShockTemperatureFloor", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ShockTemperatureFloor);
        return rv;
    }

    if (strncmp(parameter_name, "StorePreShockFields", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) StorePreShockFields);
        return rv;
    }

    if (strncmp(parameter_name, "FindShocksOnlyOnOutput", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) FindShocksOnlyOnOutput);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationFieldType", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiationFieldType);
        return rv;
    }

    if (strncmp(parameter_name, "AdjustUVBackground", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) AdjustUVBackground);
        return rv;
    }

    if (strncmp(parameter_name, "AdjustUVBackgroundHighRedshift", 31) == 0) {
        PyObject *rv = PyLong_FromLong((long) AdjustUVBackgroundHighRedshift);
        return rv;
    }

    if (strncmp(parameter_name, "SetUVBAmplitude", 16) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SetUVBAmplitude);
        return rv;
    }

    if (strncmp(parameter_name, "SetHeIIHeatingScale", 20) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SetHeIIHeatingScale);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationFieldLevelRecompute", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiationFieldLevelRecompute);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationXRaySecondaryIon", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiationXRaySecondaryIon);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationXRayComptonHeating", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiationXRayComptonHeating);
        return rv;
    }

    if (strncmp(parameter_name, "TabulatedLWBackground", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) TabulatedLWBackground);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationFieldRedshift", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiationFieldRedshift);
        return rv;
    }

    if (strncmp(parameter_name, "PhotoelectricHeating", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) PhotoelectricHeating);
        return rv;
    }

    if (strncmp(parameter_name, "PhotoelectricHeatingRate", 25) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) PhotoelectricHeatingRate);
        return rv;
    }

    if (strncmp(parameter_name, "OutputCoolingTime", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputCoolingTime);
        return rv;
    }

    if (strncmp(parameter_name, "OutputTemperature", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputTemperature);
        return rv;
    }

    if (strncmp(parameter_name, "OutputDustTemperature", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputDustTemperature);
        return rv;
    }

    if (strncmp(parameter_name, "OutputSmoothedDarkMatter", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputSmoothedDarkMatter);
        return rv;
    }

    if (strncmp(parameter_name, "SmoothedDarkMatterNeighbors", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) SmoothedDarkMatterNeighbors);
        return rv;
    }

    if (strncmp(parameter_name, "OutputGriddedStarParticle", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputGriddedStarParticle);
        return rv;
    }

    if (strncmp(parameter_name, "ZEUSLinearArtificialViscosity", 30) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ZEUSLinearArtificialViscosity);
        return rv;
    }

    if (strncmp(parameter_name, "ZEUSQuadraticArtificialViscosity", 33) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ZEUSQuadraticArtificialViscosity);
        return rv;
    }

    if (strncmp(parameter_name, "UseMinimumPressureSupport", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseMinimumPressureSupport);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumPressureSupportParameter", 32) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MinimumPressureSupportParameter);
        return rv;
    }

    if (strncmp(parameter_name, "RefineRegionTimeType", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) RefineRegionTimeType);
        return rv;
    }

    if (strncmp(parameter_name, "EvolveRefineRegionNtimes", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) EvolveRefineRegionNtimes);
        return rv;
    }

    if (strncmp(parameter_name, "MustRefineRegionTimeType", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) MustRefineRegionTimeType);
        return rv;
    }

    if (strncmp(parameter_name, "EvolveMustRefineRegionNtimes", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) EvolveMustRefineRegionNtimes);
        return rv;
    }

    if (strncmp(parameter_name, "UseCoolingRefineRegion", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseCoolingRefineRegion);
        return rv;
    }

    if (strncmp(parameter_name, "EvolveCoolingRefineRegion", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) EvolveCoolingRefineRegion);
        return rv;
    }

    if (strncmp(parameter_name, "CoolingRefineRegionTimeType", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) CoolingRefineRegionTimeType);
        return rv;
    }

    if (strncmp(parameter_name, "EvolveCoolingRefineRegionNtimes", 32) == 0) {
        PyObject *rv = PyLong_FromLong((long) EvolveCoolingRefineRegionNtimes);
        return rv;
    }

    if (strncmp(parameter_name, "MyProcessorNumber", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) MyProcessorNumber);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfProcessors", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfProcessors);
        return rv;
    }

    if (strncmp(parameter_name, "CommunicationTime", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CommunicationTime);
        return rv;
    }

    if (strncmp(parameter_name, "ParallelRootGridIO", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ParallelRootGridIO);
        return rv;
    }

    if (strncmp(parameter_name, "ParallelParticleIO", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ParallelParticleIO);
        return rv;
    }

    if (strncmp(parameter_name, "Unigrid", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) Unigrid);
        return rv;
    }

    if (strncmp(parameter_name, "CubeDumpEnabled", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) CubeDumpEnabled);
        return rv;
    }

    if (strncmp(parameter_name, "PartitionNestedGrids", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) PartitionNestedGrids);
        return rv;
    }

    if (strncmp(parameter_name, "StaticPartitionNestedGrids", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) StaticPartitionNestedGrids);
        return rv;
    }

    if (strncmp(parameter_name, "ExtractFieldsOnly", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExtractFieldsOnly);
        return rv;
    }

    if (strncmp(parameter_name, "First_Pass", 11) == 0) {
        PyObject *rv = PyLong_FromLong((long) First_Pass);
        return rv;
    }

    if (strncmp(parameter_name, "UnigridTranspose", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) UnigridTranspose);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfRootGridTilesPerDimensionPerProcessor", 46) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfRootGridTilesPerDimensionPerProcessor);
        return rv;
    }

    if (strncmp(parameter_name, "CosmologySimulationNumberOfInitialGrids", 40) == 0) {
        PyObject *rv = PyLong_FromLong((long) CosmologySimulationNumberOfInitialGrids);
        return rv;
    }

    if (strncmp(parameter_name, "OutputOnDensity", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputOnDensity);
        return rv;
    }

    if (strncmp(parameter_name, "StartDensityOutputs", 20) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StartDensityOutputs);
        return rv;
    }

    if (strncmp(parameter_name, "CurrentDensityOutput", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CurrentDensityOutput);
        return rv;
    }

    if (strncmp(parameter_name, "CurrentMaximumDensity", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CurrentMaximumDensity);
        return rv;
    }

    if (strncmp(parameter_name, "IncrementDensityOutput", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) IncrementDensityOutput);
        return rv;
    }

    if (strncmp(parameter_name, "StopFirstTimeAtDensity", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StopFirstTimeAtDensity);
        return rv;
    }

    if (strncmp(parameter_name, "StopFirstTimeAtMetalEnrichedDensity", 36) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StopFirstTimeAtMetalEnrichedDensity);
        return rv;
    }

    if (strncmp(parameter_name, "CurrentMaximumMetalEnrichedDensity", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CurrentMaximumMetalEnrichedDensity);
        return rv;
    }

    if (strncmp(parameter_name, "EnrichedMetalFraction", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) EnrichedMetalFraction);
        return rv;
    }

    if (strncmp(parameter_name, "PythonTopGridSkip", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) PythonTopGridSkip);
        return rv;
    }

    if (strncmp(parameter_name, "PythonSubcycleSkip", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) PythonSubcycleSkip);
        return rv;
    }

    if (strncmp(parameter_name, "PythonReloadScript", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) PythonReloadScript);
        return rv;
    }

    if (strncmp(parameter_name, "InlineHaloFinder", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) InlineHaloFinder);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderSubfind", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) HaloFinderSubfind);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderOutputParticleList", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) HaloFinderOutputParticleList);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderMinimumSize", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) HaloFinderMinimumSize);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderCycleSkip", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) HaloFinderCycleSkip);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderRunAfterOutput", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) HaloFinderRunAfterOutput);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderLinkingLength", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) HaloFinderLinkingLength);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderTimestep", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) HaloFinderTimestep);
        return rv;
    }

    if (strncmp(parameter_name, "HaloFinderLastTime", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) HaloFinderLastTime);
        return rv;
    }

    if (strncmp(parameter_name, "DepositPositionsParticleSmoothRadius", 37) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DepositPositionsParticleSmoothRadius);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumPressureJumpForRefinement", 33) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MinimumPressureJumpForRefinement);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumEnergyRatioForRefinement", 32) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MinimumEnergyRatioForRefinement);
        return rv;
    }

    if (strncmp(parameter_name, "RefineByJeansLengthSafetyFactor", 32) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RefineByJeansLengthSafetyFactor);
        return rv;
    }

    if (strncmp(parameter_name, "JeansRefinementColdTemperature", 31) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) JeansRefinementColdTemperature);
        return rv;
    }

    if (strncmp(parameter_name, "MustRefineParticlesRefineToLevel", 33) == 0) {
        PyObject *rv = PyLong_FromLong((long) MustRefineParticlesRefineToLevel);
        return rv;
    }

    if (strncmp(parameter_name, "MustRefineParticlesRefineToLevelAutoAdjust", 43) == 0) {
        PyObject *rv = PyLong_FromLong((long) MustRefineParticlesRefineToLevelAutoAdjust);
        return rv;
    }

    if (strncmp(parameter_name, "MustRefineParticlesMinimumMass", 31) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MustRefineParticlesMinimumMass);
        return rv;
    }

    if (strncmp(parameter_name, "MustRefineParticlesCreateParticles", 35) == 0) {
        PyObject *rv = PyLong_FromLong((long) MustRefineParticlesCreateParticles);
        return rv;
    }

    if (strncmp(parameter_name, "MinimumShearForRefinement", 26) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MinimumShearForRefinement);
        return rv;
    }

    if (strncmp(parameter_name, "OldShearMethod", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) OldShearMethod);
        return rv;
    }

    if (strncmp(parameter_name, "RefineByResistiveLengthSafetyFactor", 36) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RefineByResistiveLengthSafetyFactor);
        return rv;
    }

    if (strncmp(parameter_name, "ShockwaveRefinementMinMach", 27) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ShockwaveRefinementMinMach);
        return rv;
    }

    if (strncmp(parameter_name, "ShockwaveRefinementMinVelocity", 31) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ShockwaveRefinementMinVelocity);
        return rv;
    }

    if (strncmp(parameter_name, "ShockwaveRefinementMaxLevel", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShockwaveRefinementMaxLevel);
        return rv;
    }

    if (strncmp(parameter_name, "SecondDerivativeEpsilon", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SecondDerivativeEpsilon);
        return rv;
    }

    if (strncmp(parameter_name, "NohProblemFullBox", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) NohProblemFullBox);
        return rv;
    }

    if (strncmp(parameter_name, "ComovingCoordinates", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) ComovingCoordinates);
        return rv;
    }

    if (strncmp(parameter_name, "StarParticleCreation", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) StarParticleCreation);
        return rv;
    }

    if (strncmp(parameter_name, "StarParticleFeedback", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) StarParticleFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "StarParticleRadiativeFeedback", 30) == 0) {
        PyObject *rv = PyLong_FromLong((long) StarParticleRadiativeFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfParticleAttributes", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfParticleAttributes);
        return rv;
    }

    if (strncmp(parameter_name, "AddParticleAttributes", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) AddParticleAttributes);
        return rv;
    }

    if (strncmp(parameter_name, "BigStarFormation", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) BigStarFormation);
        return rv;
    }

    if (strncmp(parameter_name, "BigStarFormationDone", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) BigStarFormationDone);
        return rv;
    }

    if (strncmp(parameter_name, "BigStarSeparation", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) BigStarSeparation);
        return rv;
    }

    if (strncmp(parameter_name, "SimpleRampTime", 15) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SimpleRampTime);
        return rv;
    }

    if (strncmp(parameter_name, "StarFormationOncePerRootGridTimeStep", 37) == 0) {
        PyObject *rv = PyLong_FromLong((long) StarFormationOncePerRootGridTimeStep);
        return rv;
    }

    if (strncmp(parameter_name, "CubeDumps", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) CubeDumps);
        return rv;
    }

    if (strncmp(parameter_name, "TracerParticleOn", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) TracerParticleOn);
        return rv;
    }

    if (strncmp(parameter_name, "TracerParticleOutputVelocity", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) TracerParticleOutputVelocity);
        return rv;
    }

    if (strncmp(parameter_name, "TracerParticleCreationSpacing", 30) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) TracerParticleCreationSpacing);
        return rv;
    }

    if (strncmp(parameter_name, "ParticleTypeInFile", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ParticleTypeInFile);
        return rv;
    }

    if (strncmp(parameter_name, "OutputParticleTypeGrouping", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputParticleTypeGrouping);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalBoundaryIO", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExternalBoundaryIO);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalBoundaryTypeIO", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExternalBoundaryTypeIO);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalBoundaryValueIO", 24) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExternalBoundaryValueIO);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalBoundaryField", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExternalBoundaryField);
        return rv;
    }

    if (strncmp(parameter_name, "SimpleConstantBoundary", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) SimpleConstantBoundary);
        return rv;
    }

    if (strncmp(parameter_name, "StoreDomainBoundaryMassFlux", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) StoreDomainBoundaryMassFlux);
        return rv;
    }

    if (strncmp(parameter_name, "LoadGridDataAtStart", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) LoadGridDataAtStart);
        return rv;
    }

    if (strncmp(parameter_name, "PrevParameterFileName", 22) == 0) {
        PyObject *rv = PyString_FromString((char *) PrevParameterFileName);
        return rv;
    }
#ifdef MPI_INSTRUMENTATION

    if (strncmp(parameter_name, "tracename", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) tracename);
        return rv;
    }

    if (strncmp(parameter_name, "flagging_count", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) flagging_count);
        return rv;
    }

    if (strncmp(parameter_name, "in_count", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) in_count);
        return rv;
    }

    if (strncmp(parameter_name, "out_count", 10) == 0) {
        PyObject *rv = PyLong_FromLong((long) out_count);
        return rv;
    }

    if (strncmp(parameter_name, "moving_count", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) moving_count);
        return rv;
    }

    if (strncmp(parameter_name, "flagging_pct", 13) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) flagging_pct);
        return rv;
    }

    if (strncmp(parameter_name, "moving_pct", 11) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) moving_pct);
        return rv;
    }
#endif /* MPI_INSTRUMENTATION */ 

    if (strncmp(parameter_name, "name", 5) == 0) {
        PyObject *rv = PyString_FromString((char *) name);
        return rv;
    }

    if (strncmp(parameter_name, "traceMPI", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) traceMPI);
        return rv;
    }
#ifdef MEM_TRACE

    if (strncmp(parameter_name, "traceMEM", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) traceMEM);
        return rv;
    }

    if (strncmp(parameter_name, "memtracename", 13) == 0) {
        PyObject *rv = PyString_FromString((char *) memtracename);
        return rv;
    }
#endif /* MEM_TRACE */ 

    if (strncmp(parameter_name, "MovieSkipTimestep", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) MovieSkipTimestep);
        return rv;
    }

    if (strncmp(parameter_name, "Movie3DVolumes", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) Movie3DVolumes);
        return rv;
    }

    if (strncmp(parameter_name, "MovieVertexCentered", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) MovieVertexCentered);
        return rv;
    }

    if (strncmp(parameter_name, "NewMovieDumpNumber", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) NewMovieDumpNumber);
        return rv;
    }

    if (strncmp(parameter_name, "NewMovieParticleOn", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) NewMovieParticleOn);
        return rv;
    }

    if (strncmp(parameter_name, "UseHydro", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseHydro);
        return rv;
    }

    if (strncmp(parameter_name, "Coordinate", 11) == 0) {
        PyObject *rv = PyLong_FromLong((long) Coordinate);
        return rv;
    }

    if (strncmp(parameter_name, "NSpecies", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) NSpecies);
        return rv;
    }

    if (strncmp(parameter_name, "NColor", 7) == 0) {
        PyObject *rv = PyLong_FromLong((long) NColor);
        return rv;
    }

    if (strncmp(parameter_name, "Theta_Limiter", 14) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) Theta_Limiter);
        return rv;
    }

    if (strncmp(parameter_name, "UsePhysicalUnit", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) UsePhysicalUnit);
        return rv;
    }

    if (strncmp(parameter_name, "iden", 5) == 0) {
        PyObject *rv = PyLong_FromLong((long) iden);
        return rv;
    }

    if (strncmp(parameter_name, "ietot", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) ietot);
        return rv;
    }

    if (strncmp(parameter_name, "ivx", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) ivx);
        return rv;
    }

    if (strncmp(parameter_name, "ivy", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) ivy);
        return rv;
    }

    if (strncmp(parameter_name, "ivz", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) ivz);
        return rv;
    }

    if (strncmp(parameter_name, "iBx", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iBx);
        return rv;
    }

    if (strncmp(parameter_name, "iBy", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iBy);
        return rv;
    }

    if (strncmp(parameter_name, "iBz", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iBz);
        return rv;
    }

    if (strncmp(parameter_name, "iPhi", 5) == 0) {
        PyObject *rv = PyLong_FromLong((long) iPhi);
        return rv;
    }

    if (strncmp(parameter_name, "ieint", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) ieint);
        return rv;
    }

    if (strncmp(parameter_name, "iD", 3) == 0) {
        PyObject *rv = PyLong_FromLong((long) iD);
        return rv;
    }

    if (strncmp(parameter_name, "iEtot", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) iEtot);
        return rv;
    }

    if (strncmp(parameter_name, "iS1", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iS1);
        return rv;
    }

    if (strncmp(parameter_name, "iS2", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iS2);
        return rv;
    }

    if (strncmp(parameter_name, "iS3", 4) == 0) {
        PyObject *rv = PyLong_FromLong((long) iS3);
        return rv;
    }

    if (strncmp(parameter_name, "iEint", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) iEint);
        return rv;
    }

    if (strncmp(parameter_name, "SmallRho", 9) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmallRho);
        return rv;
    }

    if (strncmp(parameter_name, "SmallP", 7) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmallP);
        return rv;
    }

    if (strncmp(parameter_name, "SmallEint", 10) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmallEint);
        return rv;
    }

    if (strncmp(parameter_name, "SmallT", 7) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmallT);
        return rv;
    }

    if (strncmp(parameter_name, "MaximumAlvenSpeed", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MaximumAlvenSpeed);
        return rv;
    }

    if (strncmp(parameter_name, "NEQ_HYDRO", 10) == 0) {
        PyObject *rv = PyLong_FromLong((long) NEQ_HYDRO);
        return rv;
    }

    if (strncmp(parameter_name, "NEQ_MHD", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) NEQ_MHD);
        return rv;
    }

    if (strncmp(parameter_name, "ReconstructionMethod", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) ReconstructionMethod);
        return rv;
    }

    if (strncmp(parameter_name, "PositiveReconstruction", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) PositiveReconstruction);
        return rv;
    }

    if (strncmp(parameter_name, "RiemannSolverFallback", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) RiemannSolverFallback);
        return rv;
    }

    if (strncmp(parameter_name, "RiemannSolver", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) RiemannSolver);
        return rv;
    }

    if (strncmp(parameter_name, "ConservativeReconstruction", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) ConservativeReconstruction);
        return rv;
    }

    if (strncmp(parameter_name, "EOSType", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) EOSType);
        return rv;
    }

    if (strncmp(parameter_name, "EOSSoundSpeed", 14) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) EOSSoundSpeed);
        return rv;
    }

    if (strncmp(parameter_name, "EOSCriticalDensity", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) EOSCriticalDensity);
        return rv;
    }

    if (strncmp(parameter_name, "EOSGamma", 9) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) EOSGamma);
        return rv;
    }

    if (strncmp(parameter_name, "C_h", 4) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) C_h);
        return rv;
    }

    if (strncmp(parameter_name, "C_p", 4) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) C_p);
        return rv;
    }

    if (strncmp(parameter_name, "DivBDampingLength", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DivBDampingLength);
        return rv;
    }

    if (strncmp(parameter_name, "UseConstantAcceleration", 24) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseConstantAcceleration);
        return rv;
    }

    if (strncmp(parameter_name, "Mu", 3) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) Mu);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalGravity", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) ExternalGravity);
        return rv;
    }

    if (strncmp(parameter_name, "StringKick", 11) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StringKick);
        return rv;
    }

    if (strncmp(parameter_name, "StringKickDimension", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) StringKickDimension);
        return rv;
    }

    if (strncmp(parameter_name, "UseFloor", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseFloor);
        return rv;
    }

    if (strncmp(parameter_name, "UseViscosity", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseViscosity);
        return rv;
    }

    if (strncmp(parameter_name, "ViscosityCoefficient", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ViscosityCoefficient);
        return rv;
    }

    if (strncmp(parameter_name, "UseAmbipolarDiffusion", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseAmbipolarDiffusion);
        return rv;
    }

    if (strncmp(parameter_name, "UseResistivity", 15) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseResistivity);
        return rv;
    }

    if (strncmp(parameter_name, "HaloConcentration", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) HaloConcentration);
        return rv;
    }

    if (strncmp(parameter_name, "HaloRedshift", 13) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) HaloRedshift);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalGravityConstant", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ExternalGravityConstant);
        return rv;
    }

    if (strncmp(parameter_name, "ExternalGravityDensity", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ExternalGravityDensity);
        return rv;
    }

    if (strncmp(parameter_name, "UsePoissonDivergenceCleaning", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) UsePoissonDivergenceCleaning);
        return rv;
    }

    if (strncmp(parameter_name, "PoissonDivergenceCleaningBoundaryBuffer", 40) == 0) {
        PyObject *rv = PyLong_FromLong((long) PoissonDivergenceCleaningBoundaryBuffer);
        return rv;
    }

    if (strncmp(parameter_name, "PoissonDivergenceCleaningThreshold", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) PoissonDivergenceCleaningThreshold);
        return rv;
    }

    if (strncmp(parameter_name, "PoissonApproximationThreshold", 30) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) PoissonApproximationThreshold);
        return rv;
    }

    if (strncmp(parameter_name, "PoissonBoundaryType", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) PoissonBoundaryType);
        return rv;
    }

    if (strncmp(parameter_name, "ShiningParticleID", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShiningParticleID);
        return rv;
    }

    if (strncmp(parameter_name, "SinkMergeDistance", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SinkMergeDistance);
        return rv;
    }

    if (strncmp(parameter_name, "SinkMergeMass", 14) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SinkMergeMass);
        return rv;
    }

    if (strncmp(parameter_name, "TotalSinkMass", 14) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) TotalSinkMass);
        return rv;
    }

    if (strncmp(parameter_name, "StellarWindFeedback", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) StellarWindFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "StellarWindTurnOnMass", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StellarWindTurnOnMass);
        return rv;
    }

    if (strncmp(parameter_name, "MSStellarWindTurnOnMass", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MSStellarWindTurnOnMass);
        return rv;
    }

    if (strncmp(parameter_name, "NBodyDirectSummation", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) NBodyDirectSummation);
        return rv;
    }

    if (strncmp(parameter_name, "UseDrivingField", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseDrivingField);
        return rv;
    }

    if (strncmp(parameter_name, "DrivingEfficiency", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) DrivingEfficiency);
        return rv;
    }

    if (strncmp(parameter_name, "UseCUDA", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseCUDA);
        return rv;
    }

    if (strncmp(parameter_name, "ran1_init", 10) == 0) {
        PyObject *rv = PyLong_FromLong((long) ran1_init);
        return rv;
    }

    if (strncmp(parameter_name, "rand_init", 10) == 0) {
        PyObject *rv = PyLong_FromLong((long) rand_init);
        return rv;
    }
#ifdef STAGE_INPUT

    if (strncmp(parameter_name, "StageInput", 11) == 0) {
        PyObject *rv = PyLong_FromLong((long) StageInput);
        return rv;
    }

    if (strncmp(parameter_name, "LocalPath", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) LocalPath);
        return rv;
    }

    if (strncmp(parameter_name, "GlobalPath", 11) == 0) {
        PyObject *rv = PyString_FromString((char *) GlobalPath);
        return rv;
    }
#endif /* STAGE_INPUT */ 
#ifdef USE_PYTHON

    if (strncmp(parameter_name, "NumberOfPythonCalls", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfPythonCalls);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfPythonTopGridCalls", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfPythonTopGridCalls);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfPythonSubcycleCalls", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfPythonSubcycleCalls);
        return rv;
    }
#endif /* USE_PYTHON */ 

    if (strncmp(parameter_name, "MetalCooling", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) MetalCooling);
        return rv;
    }

    if (strncmp(parameter_name, "CIECooling", 11) == 0) {
        PyObject *rv = PyLong_FromLong((long) CIECooling);
        return rv;
    }

    if (strncmp(parameter_name, "H2OpticalDepthApproximation", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) H2OpticalDepthApproximation);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransfer", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransfer);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferHydrogenOnly", 30) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferHydrogenOnly);
        return rv;
    }
#ifdef TRANSFER

    if (strncmp(parameter_name, "PhotonTime", 11) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) PhotonTime);
        return rv;
    }

    if (strncmp(parameter_name, "dtPhoton", 9) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) dtPhoton);
        return rv;
    }
#endif /* TRANSFER */ 
#ifdef DEFINE_STORAGE
#endif /* DEFINE_STORAGE */ 

    if (strncmp(parameter_name, "RadiativeTransferSourceRadius", 30) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferSourceRadius);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferPropagationSpeedFraction", 42) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferPropagationSpeedFraction);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferPropagationDistance", 37) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferPropagationDistance);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferSplitPhotonRadius", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferSplitPhotonRadius);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferRaysPerCell", 29) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferRaysPerCell);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferInitialHEALPixLevel", 37) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferInitialHEALPixLevel);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferPhotonEscapeRadius", 36) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferPhotonEscapeRadius);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferInterpolateField", 34) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferInterpolateField);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferTimestepVelocityLimit", 39) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferTimestepVelocityLimit);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferTimestepVelocityLevel", 39) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferTimestepVelocityLevel);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferSourceClustering", 34) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferSourceClustering);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferPhotonMergeRadius", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferPhotonMergeRadius);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationPressure", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiationPressure);
        return rv;
    }

    if (strncmp(parameter_name, "RadiationPressureScale", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiationPressureScale);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferOpticallyThinH2", 33) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferOpticallyThinH2);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferOpticallyThinH2CharLength", 43) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferOpticallyThinH2CharLength);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferPeriodicBoundary", 34) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferPeriodicBoundary);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferFLDCallOnLevel", 32) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferFLDCallOnLevel);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferHIIRestrictedTimestep", 39) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferHIIRestrictedTimestep);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferAdaptiveTimestep", 34) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferAdaptiveTimestep);
        return rv;
    }

    if (strncmp(parameter_name, "GlobalMaximumkphIfront", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GlobalMaximumkphIfront);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferSourceBeamAngle", 33) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferSourceBeamAngle);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferTraceSpectrum", 31) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferTraceSpectrum);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferLoadBalance", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferLoadBalance);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferFluxBackgroundLimit", 37) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferFluxBackgroundLimit);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferRayMaximumLength", 34) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferRayMaximumLength);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferUseH2Shielding", 32) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferUseH2Shielding);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferH2ShieldType", 30) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferH2ShieldType);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferH2IIDiss", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferH2IIDiss);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferHubbleTimeFraction", 36) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) RadiativeTransferHubbleTimeFraction);
        return rv;
    }
#ifdef MEMORY_POOL
#endif /* MEMORY_POOL */ 

    if (strncmp(parameter_name, "ConductionDynamicRebuildHierarchy", 34) == 0) {
        PyObject *rv = PyLong_FromLong((long) ConductionDynamicRebuildHierarchy);
        return rv;
    }

    if (strncmp(parameter_name, "ConductionDynamicRebuildMinLevel", 33) == 0) {
        PyObject *rv = PyLong_FromLong((long) ConductionDynamicRebuildMinLevel);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferCoupledRateSolver", 35) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferCoupledRateSolver);
        return rv;
    }

    if (strncmp(parameter_name, "RadiativeTransferFLD", 21) == 0) {
        PyObject *rv = PyLong_FromLong((long) RadiativeTransferFLD);
        return rv;
    }

    if (strncmp(parameter_name, "ImplicitProblem", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) ImplicitProblem);
        return rv;
    }

    if (strncmp(parameter_name, "StarMakerEmissivityField", 25) == 0) {
        PyObject *rv = PyLong_FromLong((long) StarMakerEmissivityField);
        return rv;
    }

    if (strncmp(parameter_name, "uv_param", 9) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) uv_param);
        return rv;
    }

    if (strncmp(parameter_name, "AngularVelocity", 16) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) AngularVelocity);
        return rv;
    }

    if (strncmp(parameter_name, "VelocityGradient", 17) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) VelocityGradient);
        return rv;
    }

    if (strncmp(parameter_name, "ShearingBoundaryDirection", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShearingBoundaryDirection);
        return rv;
    }

    if (strncmp(parameter_name, "ShearingVelocityDirection", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShearingVelocityDirection);
        return rv;
    }

    if (strncmp(parameter_name, "ShearingOtherDirection", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShearingOtherDirection);
        return rv;
    }

    if (strncmp(parameter_name, "UseMHD", 7) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseMHD);
        return rv;
    }

    if (strncmp(parameter_name, "MaxVelocityIndex", 17) == 0) {
        PyObject *rv = PyLong_FromLong((long) MaxVelocityIndex);
        return rv;
    }

    if (strncmp(parameter_name, "ShearingBoxProblemType", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) ShearingBoxProblemType);
        return rv;
    }

    if (strncmp(parameter_name, "IsothermalSoundSpeed", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) IsothermalSoundSpeed);
        return rv;
    }

    if (strncmp(parameter_name, "MoveParticlesBetweenSiblings", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) MoveParticlesBetweenSiblings);
        return rv;
    }

    if (strncmp(parameter_name, "ParticleSplitterIterations", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) ParticleSplitterIterations);
        return rv;
    }

    if (strncmp(parameter_name, "ParticleSplitterChildrenParticleSeparation", 43) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ParticleSplitterChildrenParticleSeparation);
        return rv;
    }

    if (strncmp(parameter_name, "ParticleSplitterRandomSeed", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) ParticleSplitterRandomSeed);
        return rv;
    }

    if (strncmp(parameter_name, "ParticleSplitterMustRefine", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) ParticleSplitterMustRefine);
        return rv;
    }

    if (strncmp(parameter_name, "ResetMagneticField", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) ResetMagneticField);
        return rv;
    }

    if (strncmp(parameter_name, "MBHParticleIO", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) MBHParticleIO);
        return rv;
    }

    if (strncmp(parameter_name, "OutputWhenJetsHaveNotEjected", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) OutputWhenJetsHaveNotEjected);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfRadiationParticles", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfRadiationParticles);
        return rv;
    }

    if (strncmp(parameter_name, "NumberOfEnergyBins", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) NumberOfEnergyBins);
        return rv;
    }

    if (strncmp(parameter_name, "VelAnyl", 8) == 0) {
        PyObject *rv = PyLong_FromLong((long) VelAnyl);
        return rv;
    }

    if (strncmp(parameter_name, "BAnyl", 6) == 0) {
        PyObject *rv = PyLong_FromLong((long) BAnyl);
        return rv;
    }

    if (strncmp(parameter_name, "WriteExternalAccel", 19) == 0) {
        PyObject *rv = PyLong_FromLong((long) WriteExternalAccel);
        return rv;
    }

    if (strncmp(parameter_name, "UseGasDrag", 11) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseGasDrag);
        return rv;
    }

    if (strncmp(parameter_name, "GasDragCoefficient", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GasDragCoefficient);
        return rv;
    }

    if (strncmp(parameter_name, "current_error", 14) == 0) {
        PyObject *rv = PyString_FromString((char *) current_error);
        return rv;
    }

    if (strncmp(parameter_name, "IsotropicConduction", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) IsotropicConduction);
        return rv;
    }

    if (strncmp(parameter_name, "AnisotropicConduction", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) AnisotropicConduction);
        return rv;
    }

    if (strncmp(parameter_name, "IsotropicConductionSpitzerFraction", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) IsotropicConductionSpitzerFraction);
        return rv;
    }

    if (strncmp(parameter_name, "AnisotropicConductionSpitzerFraction", 37) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) AnisotropicConductionSpitzerFraction);
        return rv;
    }

    if (strncmp(parameter_name, "ConductionCourantSafetyNumber", 30) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ConductionCourantSafetyNumber);
        return rv;
    }

    if (strncmp(parameter_name, "SpeedOfLightTimeStepLimit", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) SpeedOfLightTimeStepLimit);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHFeedback", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) ClusterSMBHFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHJetMdot", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHJetMdot);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHJetVelocity", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHJetVelocity);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHJetRadius", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHJetRadius);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHJetLaunchOffset", 27) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHJetLaunchOffset);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHStartTime", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHStartTime);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHTramp", 17) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHTramp);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHJetOpenAngleRadius", 30) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHJetOpenAngleRadius);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHFastJetRadius", 25) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHFastJetRadius);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHFastJetVelocity", 27) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHFastJetVelocity);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHJetEdot", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHJetEdot);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHKineticFraction", 27) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHKineticFraction);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHJetAngleTheta", 25) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHJetAngleTheta);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHJetAnglePhi", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHJetAnglePhi);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHJetPrecessionPeriod", 31) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHJetPrecessionPeriod);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHCalculateGasMass", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) ClusterSMBHCalculateGasMass);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHFeedbackSwitch", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) ClusterSMBHFeedbackSwitch);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHEnoughColdGas", 25) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHEnoughColdGas);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHAccretionTime", 25) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHAccretionTime);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHJetDim", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) ClusterSMBHJetDim);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHAccretionEpsilon", 28) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHAccretionEpsilon);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHDiskRadius", 22) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHDiskRadius);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHBCG", 15) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHBCG);
        return rv;
    }

    if (strncmp(parameter_name, "ClusterSMBHMass", 16) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ClusterSMBHMass);
        return rv;
    }

    if (strncmp(parameter_name, "EllipticalGalaxyRe", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) EllipticalGalaxyRe);
        return rv;
    }

    if (strncmp(parameter_name, "OldStarFeedbackAlpha", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) OldStarFeedbackAlpha);
        return rv;
    }

    if (strncmp(parameter_name, "SNIaFeedbackEnergy", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SNIaFeedbackEnergy);
        return rv;
    }

    if (strncmp(parameter_name, "StellarWindRadius", 18) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StellarWindRadius);
        return rv;
    }

    if (strncmp(parameter_name, "StellarWindDensity", 19) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StellarWindDensity);
        return rv;
    }

    if (strncmp(parameter_name, "StellarWindSpeed", 17) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StellarWindSpeed);
        return rv;
    }

    if (strncmp(parameter_name, "StellarWindTemperature", 23) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) StellarWindTemperature);
        return rv;
    }

    if (strncmp(parameter_name, "MHDCTSlopeLimiter", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) MHDCTSlopeLimiter);
        return rv;
    }

    if (strncmp(parameter_name, "MHDCTDualEnergyMethod", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) MHDCTDualEnergyMethod);
        return rv;
    }

    if (strncmp(parameter_name, "MHDCTPowellSource", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) MHDCTPowellSource);
        return rv;
    }

    if (strncmp(parameter_name, "MHDCTUseSpecificEnergy", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) MHDCTUseSpecificEnergy);
        return rv;
    }

    if (strncmp(parameter_name, "WriteBoundary", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) WriteBoundary);
        return rv;
    }

    if (strncmp(parameter_name, "WriteAcceleration", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) WriteAcceleration);
        return rv;
    }

    if (strncmp(parameter_name, "TracerParticlesAddToRestart", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) TracerParticlesAddToRestart);
        return rv;
    }

    if (strncmp(parameter_name, "CT_AthenaDissipation", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) CT_AthenaDissipation);
        return rv;
    }

    if (strncmp(parameter_name, "MHD_WriteElectric", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) MHD_WriteElectric);
        return rv;
    }

    if (strncmp(parameter_name, "tiny_pressure", 14) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) tiny_pressure);
        return rv;
    }

    if (strncmp(parameter_name, "MHD_CT_Method", 14) == 0) {
        PyObject *rv = PyLong_FromLong((long) MHD_CT_Method);
        return rv;
    }

    if (strncmp(parameter_name, "MHD_ProjectB", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) MHD_ProjectB);
        return rv;
    }

    if (strncmp(parameter_name, "MHD_ProjectE", 13) == 0) {
        PyObject *rv = PyLong_FromLong((long) MHD_ProjectE);
        return rv;
    }

    if (strncmp(parameter_name, "UseMHDCT", 9) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseMHDCT);
        return rv;
    }

    if (strncmp(parameter_name, "EquationOfState", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) EquationOfState);
        return rv;
    }

    if (strncmp(parameter_name, "MHDLabel", 9) == 0) {
        PyObject *rv = PyString_FromString((char *) MHDLabel);
        return rv;
    }

    if (strncmp(parameter_name, "MHDUnits", 9) == 0) {
        PyObject *rv = PyString_FromString((char *) MHDUnits);
        return rv;
    }

    if (strncmp(parameter_name, "MHDeLabel", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) MHDeLabel);
        return rv;
    }

    if (strncmp(parameter_name, "MHDeUnits", 10) == 0) {
        PyObject *rv = PyString_FromString((char *) MHDeUnits);
        return rv;
    }

    if (strncmp(parameter_name, "CorrectParentBoundaryFlux", 26) == 0) {
        PyObject *rv = PyLong_FromLong((long) CorrectParentBoundaryFlux);
        return rv;
    }

    if (strncmp(parameter_name, "EnabledActiveParticlesCount", 28) == 0) {
        PyObject *rv = PyLong_FromLong((long) EnabledActiveParticlesCount);
        return rv;
    }

    if (strncmp(parameter_name, "ActiveParticleDensityThreshold", 31) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) ActiveParticleDensityThreshold);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarFeedback", 18) == 0) {
        PyObject *rv = PyLong_FromLong((long) SmartStarFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarEddingtonCap", 22) == 0) {
        PyObject *rv = PyLong_FromLong((long) SmartStarEddingtonCap);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarBHFeedback", 20) == 0) {
        PyObject *rv = PyLong_FromLong((long) SmartStarBHFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarBHJetFeedback", 23) == 0) {
        PyObject *rv = PyLong_FromLong((long) SmartStarBHJetFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarBHThermalFeedback", 27) == 0) {
        PyObject *rv = PyLong_FromLong((long) SmartStarBHThermalFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarBHRadiativeFeedback", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) SmartStarBHRadiativeFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarStellarRadiativeFeedback", 34) == 0) {
        PyObject *rv = PyLong_FromLong((long) SmartStarStellarRadiativeFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarFeedbackEnergyCoupling", 32) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmartStarFeedbackEnergyCoupling);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarFeedbackJetsThresholdMass", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmartStarFeedbackJetsThresholdMass);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarJetVelocity", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmartStarJetVelocity);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarSpin", 14) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmartStarSpin);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarSuperEddingtonAdjustment", 34) == 0) {
        PyObject *rv = PyLong_FromLong((long) SmartStarSuperEddingtonAdjustment);
        return rv;
    }

    if (strncmp(parameter_name, "SmartStarSMSLifetime", 21) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) SmartStarSMSLifetime);
        return rv;
    }

    if (strncmp(parameter_name, "TimingCycleSkip", 16) == 0) {
        PyObject *rv = PyLong_FromLong((long) TimingCycleSkip);
        return rv;
    }

    if (strncmp(parameter_name, "GalaxySimulationRPSWind", 24) == 0) {
        PyObject *rv = PyLong_FromLong((long) GalaxySimulationRPSWind);
        return rv;
    }

    if (strncmp(parameter_name, "GalaxySimulationRPSWindShockSpeed", 34) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GalaxySimulationRPSWindShockSpeed);
        return rv;
    }

    if (strncmp(parameter_name, "GalaxySimulationRPSWindDelay", 29) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GalaxySimulationRPSWindDelay);
        return rv;
    }

    if (strncmp(parameter_name, "GalaxySimulationRPSWindDensity", 31) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GalaxySimulationRPSWindDensity);
        return rv;
    }

    if (strncmp(parameter_name, "GalaxySimulationRPSWindTotalEnergy", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GalaxySimulationRPSWindTotalEnergy);
        return rv;
    }

    if (strncmp(parameter_name, "GalaxySimulationRPSWindPressure", 32) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GalaxySimulationRPSWindPressure);
        return rv;
    }

    if (strncmp(parameter_name, "GalaxySimulationPreWindDensity", 31) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GalaxySimulationPreWindDensity);
        return rv;
    }

    if (strncmp(parameter_name, "GalaxySimulationPreWindTotalEnergy", 35) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) GalaxySimulationPreWindTotalEnergy);
        return rv;
    }

    if (strncmp(parameter_name, "UseMagneticSupernovaFeedback", 29) == 0) {
        PyObject *rv = PyLong_FromLong((long) UseMagneticSupernovaFeedback);
        return rv;
    }

    if (strncmp(parameter_name, "MagneticSupernovaRadius", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MagneticSupernovaRadius);
        return rv;
    }

    if (strncmp(parameter_name, "MagneticSupernovaDuration", 26) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MagneticSupernovaDuration);
        return rv;
    }

    if (strncmp(parameter_name, "MagneticSupernovaEnergy", 24) == 0) {
        PyObject *rv = PyFloat_FromDouble((double) MagneticSupernovaEnergy);
        return rv;
    }
