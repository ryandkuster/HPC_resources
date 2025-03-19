# Define base directories.
BASE_DIRS := src/tools src/scripts data/docs data/db data/reference data/raw derived results

# Create the base directory structure if it doesn't exist.
.PHONY: init
init:
	@for dir in $(BASE_DIRS); do \
		mkdir -p $$dir; \
	done

NAME := $(shell date +%Y_%m_%d)_$(name)

# Function to create a new analysis directory.
.PHONY: new
new:
	@if [ -z "$(name)" ]; then \
		echo "Error: Please specify a name using 'make new name=2025_02_21_test_b'"; \
		exit 1; \
	fi
	@echo $(CURRENT_DATE)
	@mkdir -p src/$(NAME)/logs
	@mkdir -p derived/$(NAME)
	@mkdir -p results/$(NAME)
	@touch src/$(NAME)/config
	@echo 'PROJ="$${PWD%/*/*}"' >> src/$(NAME)/config
	@echo 'DATA=$${PROJ}/data' >> src/$(NAME)/config
	@echo 'SRC=$${PROJ}/src' >> src/$(NAME)/config
	@echo 'ANALYSIS=$$( basename $$PWD )' >> src/$(NAME)/config
	@echo 'DERIVED=$${PROJ}/derived/$${ANALYSIS}' >> src/$(NAME)/config
	@echo 'RESULTS=$${PROJ}/results/$${ANALYSIS}' >> src/$(NAME)/config
	@touch src/$(NAME)/00_template.sbatch
	@echo '#!/bin/bash' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH -J $(name)' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH -A acf-utk0032' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH --partition=short' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH --qos=short' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH --nodes=1' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH --cpus-per-task=40' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH --mem=100G' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH --time=0-03:00:00' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH --error=logs/job.e%J' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH --output=logs/job.o%J' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH --mail-type=END,FAIL' >> src/$(NAME)/00_template.sbatch
	@echo '#SBATCH --mail-user=rkuster@utk.edu' >> src/$(NAME)/00_template.sbatch
	@echo '' >> src/$(NAME)/00_template.sbatch
	@echo 'source config' >> src/$(NAME)/00_template.sbatch
	@echo 'echo -e "$${ANALYSIS}\n"' >> src/$(NAME)/00_template.sbatch
	@echo 'echo "project  : $${PROJ}"' >> src/$(NAME)/00_template.sbatch
	@echo 'echo "data     : $${DATA}"' >> src/$(NAME)/00_template.sbatch
	@echo 'echo "src      : $${SRC}"' >> src/$(NAME)/00_template.sbatch
	@echo 'echo "results  : $${RESULTS}"' >> src/$(NAME)/00_template.sbatch
	@echo 'echo "derived  : $${DERIVED}"' >> src/$(NAME)/00_template.sbatch
	@echo "Created directories for $(NAME)"
