function nuc = nucleation_state(numberDensity,diameter)
    %NUCLEATION_STATE Nuclei parameter, Eq. (12); diameter is cubed.
    nuc.volumeRatio = pi*numberDensity*diameter^3/6;
    nuc.alphaNuc = nuc.volumeRatio/(1+nuc.volumeRatio);
end
