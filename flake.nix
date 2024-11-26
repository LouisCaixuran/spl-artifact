{
  inputs.flake-utils.url = "github:numtide/flake-utils";
  inputs.nixpkgs.url = "nixpkgs/nixos-24.05";
  inputs.phmap-src = {
    url = "github:greg7mdp/parallel-hashmap";
    flake = false;
  };

  outputs = { self, phmap-src, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem
      (system:
        let
          pkgs = import nixpkgs {
            inherit system;
            config.allowUnfree = true;
          };
          phmap = with pkgs; stdenv.mkDerivation {
            pname = "phmap";
            version = "0.1.0";
            src = phmap-src;
            nativeBuildInputs = [ cmake ];
            cmakeFlags = [ "-DPHMAP_INSTALL=true" "-DPHMAP_BUILD_TESTS=false" "-DPHMAP_BUILD_EXAMPLES=false" ];
          };
        in
        {
          devShell = with pkgs; mkShell {
            buildInputs = [
              gnumake
              llvmPackages.bintools
              llvmPackages.openmp
              clang-tools

              clang
              bison
              boost
              flex
              gputils
              zlib
              autoconf
              automake
              phmap
              perl
              bear
              minisat
              python3
              gnum4
            ];
            MIMALLOC = "${mimalloc}";
          };
        }
      );
}
