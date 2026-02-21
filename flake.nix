{
    inputs = {
        nixpkgs.url = "github:NixOS/nixpkgs/nixpkgs-unstable";
        flake-utils.url = "github:numtide/flake-utils";
    };

    outputs = {self, ...} @inputs:
        inputs.flake-utils.lib.eachDefaultSystem (system:
            let pkgs = inputs.nixpkgs.legacyPackages.${system};
            in {
                devShells.default = pkgs.mkShell {
                    buildInputs = [
                        pkgs.python314
                        pkgs.poetry
                        pkgs.pyright
                    ];

                    shellHook = ''
                        export POETRY_VIRTUALENVS_IN_PROJECT=true
                        export POETRY_VIRTUALENVS_PATH="./.venv"

                        if [ ! -f pyproject.toml ]; then
                            echo "No pyproject.toml found. Initializing Poetry project..."
                            poetry init --no-interaction --name="ferroptosis" --python="^3.14"
                        fi

                        if [ ! -d .venv ]; then
                            echo "Creating poetry virtual environment"
                            poetry install --no-root
                        fi

                        echo "Poetry env activated"
                        echo "Python: $(python --version)"

                        exec zsh -c "source '$(pwd)/.venv/bin/activate' && exec zsh"
                    '';
                };
            }
        );
}
