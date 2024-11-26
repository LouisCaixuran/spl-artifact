FROM nixos/nix:2.23.1
RUN mkdir -p /home/spl/nix
WORKDIR /home/spl/nix
COPY flake.* .
RUN nix --extra-experimental-features nix-command --extra-experimental-features flakes develop --profile /tmp/dev
COPY sdcc-cfg /home/spl/sdcc-cfg
COPY start.sh /home/spl/sdcc-cfg
COPY run.sh /home/spl/sdcc-cfg
COPY spl_g.py /home/spl/
WORKDIR /home/spl/sdcc-cfg

