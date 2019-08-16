FROM julia:1.1.1-buster

LABEL description="Julia Image to run HMC.jl" 

WORKDIR /app
ENV JULIA_PROJECT=/app
COPY *.toml /app/
COPY build_deps.jl /app/
RUN julia -O3 build_deps.jl

COPY build_hmc.jl /app/
COPY src/Hmc.jl /app/src/Hmc.jl
RUN julia -O3 build_hmc.jl

CMD ["/bin/bash"]