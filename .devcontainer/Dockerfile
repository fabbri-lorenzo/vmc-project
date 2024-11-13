FROM debian:trixie-slim

ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

COPY library-scripts/*.sh /library-scripts/

RUN \
    apt-get update && \
    apt-get -y upgrade && \
    chmod +x /library-scripts/*.sh && \
    library-scripts/provide-dependencies-ci.sh && \
    library-scripts/provide-dependencies-local.sh && \
    library-scripts/provide-user.sh ${USERNAME} ${USER_UID} ${USER_GID} && \
    rm -rf /library-scripts