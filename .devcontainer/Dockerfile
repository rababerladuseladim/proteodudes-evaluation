# syntax=docker/dockerfile:1
FROM mcr.microsoft.com/devcontainers/miniconda:3

ARG USERNAME=vscode

# Pre-create the Claude config dir and chown it to the non-root user
# *before* it becomes a named-volume mount target in devcontainer.json,
# so the volume inherits non-root ownership on first creation
# (mirrors anthropics/claude-code/.devcontainer/Dockerfile).
RUN mkdir -p /home/${USERNAME}/.claude \
    && chown -R ${USERNAME}:${USERNAME} /home/${USERNAME}/.claude

# Bake the locked "snakemake" env into the image at build time, using
# the same explicit lockfile CI installs from (conda-incubator/setup-miniconda
# + spec-file.txt in .github/workflows/main.yaml), instead of re-solving
# env.yaml on every container start.
COPY spec-file.txt /tmp/spec-file.txt
RUN conda create --name snakemake --file /tmp/spec-file.txt \
    && rm /tmp/spec-file.txt

USER ${USERNAME}
