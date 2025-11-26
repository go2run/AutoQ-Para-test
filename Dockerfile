FROM fedora:latest
RUN dnf install -y cmake make g++ gmp-devel || sudo apt update && sudo apt install -y cmake make g++ libgmp-dev
RUN mkdir -p /workspace/build
COPY ./src ./workspace
WORKDIR /workspace/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j6 experiments
ENTRYPOINT ./experiments
