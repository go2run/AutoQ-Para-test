FROM --platform=linux/amd64 fedora:latest

# 建置時可用的參數（必要時可從 devcontainer.json 覆寫）
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=1000

# 1) 安裝必要套件（不包含 gcovr / subversion / ca-certificates / lsd）
# 2) 建立非 root 使用者，並允許無密碼 sudo（方便 VS Code 運作）
RUN dnf -y update && \
    dnf install -y \
      cmake \
      make \
      gcc-c++ \
      gmp-devel \
      git git-lfs \
      pkgconfig \
      zsh \
      sudo \
      which \
      openssh-clients \
      procps \
      catch2-devel && \
    dnf clean all

    
RUN groupadd -g ${USER_GID} ${USERNAME} || true && \
    useradd -m -u ${USER_UID} -g ${USER_GID} -s /bin/zsh ${USERNAME} || true && \
    echo "${USERNAME} ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/${USERNAME} && \
    chmod 0440 /etc/sudoers.d/${USERNAME}

# 建置目錄與原始碼複製
RUN mkdir -p /workspace
COPY ./src /workspace
RUN chown -R ${USERNAME}:${USERNAME} /workspace

# 切換到非 root 使用者進行實際編譯
USER ${USERNAME}
ENV HOME=/home/${USERNAME}
WORKDIR /workspace/build

# 使用 zsh 為預設 shell 執行後續 RUN 指令（可使開發體驗一致）
SHELL ["/bin/zsh", "-c"]

# CMake 配置與編譯（以非 root 使用者執行）
RUN mkdir -p /workspace/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j$(nproc) experiments

# 開發容器建議：不要覆寫 entrypoint，讓 VS Code 遠端容器管理啟動流程；
# 若您希望容器啟動即執行可執行檔，可改用 CMD 或 ENTRYPOINT。
CMD ["/bin/zsh"]