FROM mcr.microsoft.com/devcontainers/base:jammy

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Shared with CI. Pandoc comes from Quarto, not jammy's outdated apt package.
COPY .github/scripts/system-deps.txt /tmp/system-deps.txt
RUN apt-get update \
    && xargs -a /tmp/system-deps.txt apt-get install -y --no-install-recommends \
    && apt-get install -y --no-install-recommends curl ca-certificates jq gdebi-core \
    && curl -fsSL https://rig.r-pkg.org/deb/rig.gpg \
       -o /usr/share/keyrings/rig.gpg \
    && echo "deb [signed-by=/usr/share/keyrings/rig.gpg] https://rig.r-pkg.org/deb rig main" \
       > /etc/apt/sources.list.d/rig.list \
    && apt-get update \
    && apt-get install -y --no-install-recommends r-rig \
    && rm -rf /var/lib/apt/lists/* /tmp/system-deps.txt

# A changed lockfile invalidates this layer, including a changed R.Version.
COPY renv.lock /tmp/opal/renv.lock
RUN r_version="$(jq -er '.R.Version' /tmp/opal/renv.lock)" \
    && rig add "$r_version" --without-pak --without-sysreqs \
    && rig default "$r_version" \
    && Rscript --vanilla -e 'cat( \
       "\noptions(repos = c(CRAN = \"https://packagemanager.posit.co/cran/__linux__/jammy/latest\"))\n", \
       file = file.path(R.home("etc"), "Rprofile.site"), append = TRUE)' \
    && mkdir -p /home/vscode/.local/share/renv/cache \
    && chown -R vscode:vscode /home/vscode/.local

ENV RENV_CONFIG_REPOS_OVERRIDE=https://packagemanager.posit.co/cran/__linux__/jammy/latest
