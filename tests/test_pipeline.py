"""Tests for pipeline orchestration components."""

import json
from datetime import datetime
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from seednap.config.models import PipelineConfig
from seednap.pipeline.orchestrator import PipelineOrchestrator
from seednap.pipeline.state import PipelineState, StepState, StepStatus


# ==============================================================================
# Test StepState
# ==============================================================================


class TestStepState:
    """Test StepState class."""

    def test_create_step_state(self):
        """Test creating a step state."""
        step = StepState(name="trim")

        assert step.name == "trim"
        assert step.status == StepStatus.PENDING
        assert step.started_at is None
        assert step.completed_at is None
        assert step.duration_seconds is None
        assert step.error_message is None
        assert step.outputs == {}
        assert step.metadata == {}

    def test_start_step(self):
        """Test starting a step."""
        step = StepState(name="trim")
        step.start()

        assert step.status == StepStatus.RUNNING
        assert step.started_at is not None
        assert isinstance(step.started_at, datetime)

    def test_complete_step(self):
        """Test completing a step."""
        step = StepState(name="trim")
        step.start()

        outputs = {"trimmed_dir": Path("/tmp/trimmed")}
        step.complete(outputs)

        assert step.status == StepStatus.COMPLETED
        assert step.completed_at is not None
        assert step.duration_seconds is not None
        assert step.duration_seconds >= 0
        assert step.outputs == outputs

    def test_complete_step_without_outputs(self):
        """Test completing a step without outputs."""
        step = StepState(name="trim")
        step.start()
        step.complete()

        assert step.status == StepStatus.COMPLETED
        assert step.outputs == {}

    def test_fail_step(self):
        """Test failing a step."""
        step = StepState(name="trim")
        step.start()

        error = ValueError("Something went wrong")
        step.fail(error)

        assert step.status == StepStatus.FAILED
        assert step.completed_at is not None
        assert step.duration_seconds is not None
        assert step.error_message == "Something went wrong"

    def test_skip_step(self):
        """Test skipping a step."""
        step = StepState(name="trim")
        step.skip(reason="Already trimmed")

        assert step.status == StepStatus.SKIPPED
        assert step.metadata["skip_reason"] == "Already trimmed"

    def test_skip_step_no_reason(self):
        """Test skipping a step without reason."""
        step = StepState(name="trim")
        step.skip()

        assert step.status == StepStatus.SKIPPED
        assert "skip_reason" not in step.metadata


# ==============================================================================
# Test PipelineState
# ==============================================================================


class TestPipelineState:
    """Test PipelineState class."""

    def test_create_pipeline_state(self):
        """Test creating a pipeline state."""
        state = PipelineState(marker="teleo")

        assert state.marker == "teleo"
        assert state.started_at is not None
        assert state.completed_at is None
        assert state.steps == {}
        assert state.current_step is None

    def test_add_step(self):
        """Test adding a step."""
        state = PipelineState(marker="teleo")
        step = state.add_step("trim")

        assert "trim" in state.steps
        assert step.name == "trim"
        assert step.status == StepStatus.PENDING

    def test_add_duplicate_step(self):
        """Test adding a duplicate step returns existing step."""
        state = PipelineState(marker="teleo")
        step1 = state.add_step("trim")
        step2 = state.add_step("trim")

        assert step1 is step2
        assert len(state.steps) == 1

    def test_get_step(self):
        """Test getting a step."""
        state = PipelineState(marker="teleo")
        state.add_step("trim")

        step = state.get_step("trim")
        assert step is not None
        assert step.name == "trim"

    def test_get_nonexistent_step(self):
        """Test getting a nonexistent step."""
        state = PipelineState(marker="teleo")
        step = state.get_step("nonexistent")

        assert step is None

    def test_start_step(self):
        """Test starting a step."""
        state = PipelineState(marker="teleo")
        step = state.start_step("trim")

        assert step.status == StepStatus.RUNNING
        assert state.current_step == "trim"

    def test_start_nonexistent_step_creates_it(self):
        """Test starting a nonexistent step creates it."""
        state = PipelineState(marker="teleo")
        step = state.start_step("trim")

        assert "trim" in state.steps
        assert step.status == StepStatus.RUNNING

    def test_complete_step(self):
        """Test completing a step."""
        state = PipelineState(marker="teleo")
        state.start_step("trim")

        outputs = {"trimmed_dir": Path("/tmp/trimmed")}
        state.complete_step("trim", outputs)

        step = state.get_step("trim")
        assert step.status == StepStatus.COMPLETED
        assert step.outputs == outputs
        assert state.current_step is None

    def test_fail_step(self):
        """Test failing a step."""
        state = PipelineState(marker="teleo")
        state.start_step("trim")

        error = ValueError("Trimming failed")
        state.fail_step("trim", error)

        step = state.get_step("trim")
        assert step.status == StepStatus.FAILED
        assert step.error_message == "Trimming failed"
        assert state.current_step is None

    def test_skip_step(self):
        """Test skipping a step."""
        state = PipelineState(marker="teleo")
        state.skip_step("trim", reason="Already trimmed")

        step = state.get_step("trim")
        assert step.status == StepStatus.SKIPPED

    def test_is_step_completed(self):
        """Test checking if step is completed."""
        state = PipelineState(marker="teleo")
        state.start_step("trim")
        state.complete_step("trim")

        assert state.is_step_completed("trim") is True
        assert state.is_step_completed("dada2") is False

    def test_is_step_failed(self):
        """Test checking if step failed."""
        state = PipelineState(marker="teleo")
        state.start_step("trim")
        state.fail_step("trim", "Error")

        assert state.is_step_failed("trim") is True
        assert state.is_step_failed("dada2") is False

    def test_get_completed_steps(self):
        """Test getting completed steps."""
        state = PipelineState(marker="teleo")
        state.start_step("trim")
        state.complete_step("trim")
        state.start_step("dada2")
        state.complete_step("dada2")

        completed = state.get_completed_steps()
        assert set(completed) == {"trim", "dada2"}

    def test_get_failed_steps(self):
        """Test getting failed steps."""
        state = PipelineState(marker="teleo")
        state.start_step("trim")
        state.fail_step("trim", "Error")
        state.start_step("dada2")
        state.complete_step("dada2")

        failed = state.get_failed_steps()
        assert failed == ["trim"]

    def test_get_pending_steps(self):
        """Test getting pending steps."""
        state = PipelineState(marker="teleo")
        state.add_step("trim")
        state.add_step("dada2")
        state.add_step("taxonomy")

        state.start_step("trim")
        state.complete_step("trim")

        all_steps = ["trim", "dada2", "taxonomy", "export"]
        pending = state.get_pending_steps(all_steps)

        assert set(pending) == {"dada2", "taxonomy", "export"}

    def test_can_resume(self):
        """Test checking if pipeline can be resumed."""
        state = PipelineState(marker="teleo")

        # No completed steps - cannot resume
        assert state.can_resume() is False

        # Has completed step - can resume
        state.start_step("trim")
        state.complete_step("trim")
        assert state.can_resume() is True

        # Has running step - cannot resume
        state.start_step("dada2")
        assert state.can_resume() is False

    def test_complete_pipeline(self):
        """Test completing the pipeline."""
        state = PipelineState(marker="teleo")
        state.complete_pipeline()

        assert state.completed_at is not None
        assert state.current_step is None

    def test_get_summary(self):
        """Test getting pipeline summary."""
        state = PipelineState(marker="teleo")
        state.start_step("trim")
        state.complete_step("trim")
        state.start_step("dada2")
        state.fail_step("dada2", "Error")
        state.skip_step("taxonomy")
        state.add_step("export")

        state.complete_pipeline()

        summary = state.get_summary()

        assert summary["marker"] == "teleo"
        assert summary["total_steps"] == 4
        assert summary["completed"] == 1
        assert summary["failed"] == 1
        assert summary["skipped"] == 1
        assert summary["pending"] == 1
        assert summary["total_duration_seconds"] is not None

    def test_save_and_load_state(self, tmp_path: Path):
        """Test saving and loading pipeline state."""
        state_file = tmp_path / "state.json"

        # Create and save state
        state = PipelineState(marker="teleo", config_path=Path("/tmp/config.yaml"))
        state.start_step("trim")
        state.complete_step("trim", {"trimmed_dir": Path("/tmp/trimmed")})

        state.save(state_file)

        # Load state
        loaded_state = PipelineState.load(state_file)

        assert loaded_state.marker == "teleo"
        assert "trim" in loaded_state.steps
        assert loaded_state.is_step_completed("trim")

    def test_load_nonexistent_state(self):
        """Test loading nonexistent state file."""
        with pytest.raises(FileNotFoundError):
            PipelineState.load(Path("/nonexistent/state.json"))

    def test_from_config(self):
        """Test creating state from config."""
        state = PipelineState.from_config(marker="teleo", config_path=Path("/tmp/config.yaml"))

        assert state.marker == "teleo"
        assert state.config_path == Path("/tmp/config.yaml")


# ==============================================================================
# Test PipelineOrchestrator
# ==============================================================================


class TestPipelineOrchestrator:
    """Test PipelineOrchestrator class."""

    @pytest.fixture
    def minimal_config(self, tmp_path: Path) -> PipelineConfig:
        """Create minimal pipeline config for testing."""
        config = PipelineConfig(
            version="0.1.0",
            marker={
                "name": "teleo",
                "description": "Test marker",
                "primers": {
                    "forward": "ACACCGCCCGTCACTCT",
                    "reverse": "CTTCCGGTACACTTACCATG",
                    "name": "Teleo",
                    "target": "12S rRNA",
                    "amplicon_length": [100, 200],
                },
            },
            paths={
                "raw_data": tmp_path / "raw",
                "output": tmp_path / "outputs",
                "logs": tmp_path / "logs",
            },
            taxonomy={
                "method": "dada2",
                "databases": {
                    "dada2": {
                        "all": tmp_path / "db" / "all.fasta",
                        "species": tmp_path / "db" / "species.fasta",
                    }
                },
            },
        )

        # Create necessary directories
        (tmp_path / "raw").mkdir(parents=True, exist_ok=True)
        (tmp_path / "db").mkdir(parents=True, exist_ok=True)

        return config

    def test_create_orchestrator(self, minimal_config: PipelineConfig):
        """Test creating pipeline orchestrator."""
        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)

            assert orchestrator.config.marker.name == "teleo"
            assert orchestrator.state.marker == "teleo"
            assert len(orchestrator.state.steps) == len(minimal_config.pipeline.steps)

    def test_create_orchestrator_from_yaml(self, tmp_path: Path, minimal_config: PipelineConfig):
        """Test creating orchestrator from YAML config."""
        # Save config to YAML
        config_file = tmp_path / "config.yaml"
        import yaml

        # Convert Path objects to strings and tuples to lists for YAML serialization
        def convert_paths(obj):
            if isinstance(obj, Path):
                return str(obj)
            elif isinstance(obj, dict):
                return {k: convert_paths(v) for k, v in obj.items()}
            elif isinstance(obj, (list, tuple)):
                return [convert_paths(item) for item in obj]
            return obj

        config_dict = convert_paths(minimal_config.model_dump(mode="python"))

        with open(config_file, "w") as f:
            yaml.dump(config_dict, f)

        with patch("seednap.pipeline.orchestrator.setup_logging"):
            with patch("seednap.config.loader.load_config", return_value=minimal_config):
                orchestrator = PipelineOrchestrator(config=config_file)

                assert orchestrator.config_path == config_file

    def test_create_orchestrator_resume(self, tmp_path: Path, minimal_config: PipelineConfig):
        """Test creating orchestrator with resume."""
        state_file = tmp_path / "state.json"

        # Create initial state
        state = PipelineState.from_config(marker="teleo")
        state.start_step("trim")
        state.complete_step("trim")
        state.save(state_file)

        # Resume
        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(
                config=minimal_config, state_file=state_file, resume=True
            )

            assert orchestrator.state.is_step_completed("trim")

    def test_resume_without_state_file_fails(self, minimal_config: PipelineConfig):
        """Test that resuming without state file raises error."""
        with pytest.raises(ValueError, match="Cannot resume"):
            with patch("seednap.pipeline.orchestrator.setup_logging"):
                PipelineOrchestrator(
                    config=minimal_config,
                    state_file=Path("/nonexistent/state.json"),
                    resume=True,
                )

    def test_should_run_step_skipped(self, minimal_config: PipelineConfig):
        """Test that skipped steps should not run."""
        minimal_config.pipeline.skip = ["trim"]

        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)

            assert orchestrator._should_run_step("trim") is False

    def test_should_run_step_completed(self, minimal_config: PipelineConfig):
        """Test that completed steps should not run."""
        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)
            orchestrator.state.start_step("trim")
            orchestrator.state.complete_step("trim")

            assert orchestrator._should_run_step("trim") is False

    def test_should_run_step_failed(self, minimal_config: PipelineConfig):
        """Test that failed steps should retry."""
        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)
            orchestrator.state.start_step("trim")
            orchestrator.state.fail_step("trim", "Error")

            assert orchestrator._should_run_step("trim") is True

    @patch("seednap.pipeline.orchestrator.StandardTrimmer")
    def test_run_trim_step(self, mock_trimmer_class: MagicMock, minimal_config: PipelineConfig, tmp_path: Path):
        """Test running trim step."""
        # Create sample files
        raw_dir = minimal_config.paths.raw_data
        (raw_dir / "sample1_R1.fastq.gz").write_text("@read1\nATCG\n+\nIIII\n")
        (raw_dir / "sample1_R2.fastq.gz").write_text("@read1\nGCTA\n+\nIIII\n")

        # Mock trimmer
        mock_trimmer = MagicMock()
        mock_trimmer.trim_sample.return_value = {
            "r1_output": tmp_path / "sample1.R1.fastq",
            "r2_output": tmp_path / "sample1.R2.fastq",
        }
        mock_trimmer_class.return_value = mock_trimmer

        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)
            outputs = orchestrator.run_trim()

            assert "trimmed_dir" in outputs
            assert "samples" in outputs
            assert orchestrator.state.is_step_completed("trim")
            mock_trimmer.trim_sample.assert_called_once()

    @patch("seednap.pipeline.orchestrator.Dada2Processor")
    def test_run_dada2_step(self, mock_processor_class: MagicMock, minimal_config: PipelineConfig):
        """Test running DADA2 step."""
        # Mock processor
        mock_processor = MagicMock()
        mock_processor.process.return_value = {
            "seqtab": Path("/tmp/seqtab.rds"),
            "query_fasta": Path("/tmp/query.fasta"),
            "seqtab_clean_t": Path("/tmp/seqtab_clean_t.csv"),
        }
        mock_processor_class.return_value = mock_processor

        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)
            outputs = orchestrator.run_dada2()

            assert "seqtab" in outputs
            assert orchestrator.state.is_step_completed("dada2")
            mock_processor.process.assert_called_once()

    @patch("seednap.pipeline.orchestrator.TaxonomicAssigner")
    def test_run_taxonomy_step(
        self, mock_assigner_class: MagicMock, minimal_config: PipelineConfig, tmp_path: Path
    ):
        """Test running taxonomy step."""
        # Setup DADA2 outputs in state
        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)

            # Mark DADA2 as completed with outputs
            orchestrator.state.start_step("dada2")
            orchestrator.state.complete_step(
                "dada2",
                {
                    "query_fasta": tmp_path / "query.fasta",
                    "seqtab_clean_t": tmp_path / "seqtab.csv",
                },
            )

            # Mock assigner
            mock_assigner = MagicMock()
            mock_assigner.assign_taxonomy.return_value = {
                "taxonomy": tmp_path / "taxonomy.csv",
                "final_table": tmp_path / "complete.csv",
            }
            mock_assigner_class.return_value = mock_assigner

            outputs = orchestrator.run_taxonomy()

            assert "final_table" in outputs
            assert orchestrator.state.is_step_completed("taxonomy")
            mock_assigner.assign_taxonomy.assert_called_once()

    def test_run_taxonomy_without_dada2_fails(self, minimal_config: PipelineConfig):
        """Test that taxonomy fails if DADA2 not completed."""
        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)

            with pytest.raises(ValueError, match="DADA2 step must be completed"):
                orchestrator.run_taxonomy()

    @patch("seednap.pipeline.orchestrator.GBIFFormatter")
    def test_run_export_step(
        self, mock_formatter_class: MagicMock, minimal_config: PipelineConfig, tmp_path: Path
    ):
        """Test running export step."""
        # Setup taxonomy outputs in state
        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)

            # Mark taxonomy as completed
            taxonomy_csv = tmp_path / "taxonomy.csv"
            taxonomy_csv.write_text("sequence,kingdom\nATCG,Animalia\n")

            orchestrator.state.start_step("taxonomy")
            orchestrator.state.complete_step("taxonomy", {"final_table": taxonomy_csv})

            # Mock formatter
            mock_formatter = MagicMock()
            import pandas as pd

            mock_formatter.from_dada2_rdp.return_value = pd.DataFrame(
                {"sequence": ["ATCG"], "kingdom": ["Animalia"]}
            )
            mock_formatter_class.return_value = mock_formatter

            outputs = orchestrator.run_export()

            assert "gbif_csv" in outputs
            assert orchestrator.state.is_step_completed("export")

    def test_get_sample_list(self, minimal_config: PipelineConfig):
        """Test getting sample list from raw data directory."""
        raw_dir = minimal_config.paths.raw_data

        # Create sample files
        (raw_dir / "sample1_R1.fastq.gz").touch()
        (raw_dir / "sample1_R2.fastq.gz").touch()
        (raw_dir / "sample2_R1.fastq.gz").touch()
        (raw_dir / "sample2_R2.fastq.gz").touch()

        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)
            samples = orchestrator._get_sample_list()

            assert set(samples) == {"sample1", "sample2"}

    def test_find_read_file(self, minimal_config: PipelineConfig):
        """Test finding read files."""
        raw_dir = minimal_config.paths.raw_data
        (raw_dir / "sample1_R1.fastq.gz").touch()

        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)
            read_file = orchestrator._find_read_file("sample1", "R1")

            assert read_file == raw_dir / "sample1_R1.fastq.gz"

    def test_find_read_file_not_found(self, minimal_config: PipelineConfig):
        """Test finding nonexistent read file."""
        with patch("seednap.pipeline.orchestrator.setup_logging"):
            orchestrator = PipelineOrchestrator(config=minimal_config)

            with pytest.raises(FileNotFoundError, match="Could not find R1 file"):
                orchestrator._find_read_file("nonexistent", "R1")

    def test_state_saved_after_step(self, minimal_config: PipelineConfig, tmp_path: Path):
        """Test that state is saved after each step."""
        state_file = tmp_path / "state.json"

        with patch("seednap.pipeline.orchestrator.setup_logging"):
            with patch("seednap.pipeline.orchestrator.Dada2Processor"):
                orchestrator = PipelineOrchestrator(config=minimal_config, state_file=state_file)
                orchestrator.run_dada2()

                # Check state file exists and contains the step
                assert state_file.exists()
                with open(state_file) as f:
                    state_data = json.load(f)
                    assert "dada2" in state_data["steps"]
                    assert state_data["steps"]["dada2"]["status"] == "completed"
