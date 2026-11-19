# Minimal x86-64 Linux runtime for the Linux build of ORCA.
# ORCA itself is NOT copied into the image (licence); it is bind-mounted at run time.
FROM --platform=linux/amd64 debian:bookworm-slim
RUN apt-get update \
 && apt-get install -y --no-install-recommends openmpi-bin libopenmpi3 openssh-client ca-certificates \
 && rm -rf /var/lib/apt/lists/*
ENV OMPI_MCA_plm=isolated \
    OMPI_MCA_btl_vader_single_copy_mechanism=none \
    OMPI_MCA_rmaps_base_oversubscribe=yes \
    HOME=/tmp
