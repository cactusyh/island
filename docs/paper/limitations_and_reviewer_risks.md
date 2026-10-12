# Limitations and reviewer risks

- The current J23 execution is OpenMM Reference and is intended for short
  software checks. It does not establish production scientific MD.
- The pinned LAMMPS executable lacks KSPACE; no Ewald/PPPM runtime comparison is
  claimed.
- Triclinic cells and NPT annealing are unsupported. RDF currently requires
  orthorhombic boxes.
- PCFF `.frc`, OPLS-AA, GAFF and GAFF2 support is profile-bounded; universal
  coverage and complete msi2lmp parity are not claimed.
- Physical density convergence, replica agreement, Tg and mechanical validation
  require external trajectories and protocols that are not present here.
- Full checkpoint bytes are retained only at accepted restart boundaries by the
  opt-in store. The historical J23 short-workflow schema is preserved unchanged.
- No external property-prediction or AI benchmark is part of this paper.
