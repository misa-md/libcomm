<a name="unreleased"></a>
## [Unreleased]


<a name="v0.2.0"></a>
## [v0.2.0] - 2019-07-09
### Feat
- **domain:** add feature of setting ghost size in domain builder.

### Fix
- **$compile:** change type u_int64_t to uint64_t to fix compile error.
- **domain:** move members of doubled x coordinate for class comm::Domain to class comm::BccDomain.

### Refactor
- **domain:** move comm::Domain::Builder out of class comm::Domain.

### Test
- **domain:** correct domain tests effected by commit "change ghost size to (cut_lattice + 1)".

### BREAKING CHANGE

change domain param in comm ::fwCommLocalRegion from comm::Domain to comm::BccDomain.


<a name="v0.1.0"></a>
## [v0.1.0] - 2019-05-31
### Fix
- **domain:** change ghost size to (cut_lattice + 1) instead of cut_lattice.
- **forward-comm:** use db_lattice_size_ghost instead of cut_lattice to compute region in forward communication.


<a name="v0.1.0-beta"></a>
## [v0.1.0-beta] - 2019-04-21
### Feat
- **forward-comm:** add feature of mearsured region calculating for forwarding communication.
- **forward-comm:** add feature of region calculating for forwarding communication.

### Fix
- **domain:** fix range issue in Region::isIn function.


<a name="v0.1.0-alpha"></a>
## v0.1.0-alpha - 2019-04-14
### Docs
- **readme:** add document for lib dependencies building in README.md.

### Feat
- **communication:** add feature of reversed order of dimension loop (z,y,z).
- **communication:** add generic communication function neiSendReceive and data Packer implementation.
- **domain:** add localBuild for domain.

### Fix
- **domain:** fix compile error of domain.cpp/.h, and add comm namespace.
- **domain:** double the lattice size and lattice coord in x direction.
- **pkg:** correct the mistake spell of "cmake_build" to "cmake_lib" in file pkg.yaml.

### Refactor
- **domain:** add double-x lattice size and coord member for class Domain.
- **domain:** refator doamin implementation to use lattice priority strategy (rather than measured length priority).
- **domain:** extract domain boundary to class Region.

### Test
- **domain:** add tests for domain decomposition.


[Unreleased]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.2.0...HEAD
[v0.2.0]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.1.0...v0.2.0
[v0.1.0]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.1.0-beta...v0.1.0
[v0.1.0-beta]: https://git.hpcer.dev/HPCer/CrystalMD/CrystalMD/compare/v0.1.0-alpha...v0.1.0-beta
