"""Tests for the scnado CLI."""
import re
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from scnado.cli import app


runner = CliRunner()


class TestFindBarcodesCommand:
    """Tests for the find-barcodes command."""

    def test_find_barcodes_help(self):
        """Test that help works for find-barcodes command."""
        result = runner.invoke(app, ["find-barcodes", "--help"])
        assert result.exit_code == 0
        assert "Find barcodes" in result.stdout
        assert "--r1" in result.stdout
        assert "--r2" in result.stdout
        assert "--barcodes" in result.stdout
        assert "--output-prefix" in result.stdout

    def test_find_barcodes_missing_required_args(self):
        """Test that missing required arguments causes error."""
        result = runner.invoke(app, ["find-barcodes", "--r1", "file.fastq.gz"])
        assert result.exit_code != 0

    def test_find_barcodes_without_rust_module(self):
        """Test handling of missing Rust module (normal case without mock)."""
        result = runner.invoke(
            app,
            [
                "find-barcodes",
                "--r1", "test_R1.fastq.gz",
                "--r2", "test_R2.fastq.gz",
                "--barcodes", "barcodes.csv",
                "--output-prefix", "output/sample",
            ],
        )
        
        # Should fail with ImportError if Rust module not installed (normal case)
        assert result.exit_code == 1


class TestFragmentsCommand:
    """Tests for the fragments command."""

    def test_fragments_help(self):
        """Test that help works for fragments command."""
        result = runner.invoke(app, ["fragments", "--help"])
        assert result.exit_code == 0
        assert "Extract fragments" in result.stdout
        assert "--bam" in result.stdout
        assert "--output" in result.stdout
        assert "--barcode-regex" in result.stdout

    def test_fragments_invalid_regex(self):
        """Test that invalid regex is caught."""
        result = runner.invoke(
            app,
            [
                "fragments",
                "--bam", "aligned.bam",
                "--output", "fragments.tsv.gz",
                "--barcode-regex", "(?P<invalid",  # Invalid regex
            ],
        )
        
        assert result.exit_code == 1
        assert "Invalid regex" in result.stdout

    def test_fragments_without_rust_module(self):
        """Test handling of missing Rust module (normal case without mock)."""
        result = runner.invoke(
            app,
            [
                "fragments",
                "--bam", "aligned.bam",
                "--output", "fragments.tsv.gz",
            ],
        )
        
        # Should fail with ImportError if Rust module not installed (normal case)
        assert result.exit_code == 1


class TestWorkflowInit:
    """Tests for the workflow init command."""

    def test_workflow_init_help(self):
        """Test that help works for workflow init command."""
        result = runner.invoke(app, ["workflow", "init", "--help"])
        assert result.exit_code == 0
        assert "Initialize" in result.stdout or "initialize" in result.stdout
        assert "--outdir" in result.stdout

    def test_workflow_init_creates_files(self):
        """Test that workflow init creates config files."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(
                app,
                ["workflow", "init", "--outdir", tmpdir],
            )
            
            assert result.exit_code == 0
            assert "✓ Created config template" in result.stdout
            assert "✓ Created samples template" in result.stdout
            
            # Check files exist
            config_file = Path(tmpdir) / "config" / "config.yaml"
            samples_file = Path(tmpdir) / "config" / "samples.tsv"
            
            assert config_file.exists()
            assert samples_file.exists()

    def test_workflow_init_config_content(self):
        """Test that config.yaml has correct content."""
        with tempfile.TemporaryDirectory() as tmpdir:
            runner.invoke(
                app,
                ["workflow", "init", "--outdir", tmpdir],
            )
            
            config_file = Path(tmpdir) / "config" / "config.yaml"
            content = config_file.read_text()
            
            assert "samples:" in content
            assert "enable_cat:" in content
            assert "enable_rna:" in content
            assert "bowtie2_index:" in content
            assert "scnado_container:" in content

    def test_workflow_init_samples_content(self):
        """Test that samples.tsv has correct content."""
        with tempfile.TemporaryDirectory() as tmpdir:
            runner.invoke(
                app,
                ["workflow", "init", "--outdir", tmpdir],
            )
            
            samples_file = Path(tmpdir) / "config" / "samples.tsv"
            content = samples_file.read_text()
            
            assert "sample" in content
            assert "sublibrary" in content
            assert "fastq_r1" in content
            assert "fastq_r2" in content
            assert "modality" in content

    def test_workflow_init_default_outdir(self):
        """Test that default outdir is used when not specified."""
        with tempfile.TemporaryDirectory() as tmpdir:
            import os
            original_cwd = os.getcwd()
            try:
                os.chdir(tmpdir)
                result = runner.invoke(app, ["workflow", "init"])
                
                assert result.exit_code == 0
                project_dir = Path(tmpdir) / "scnado_project"
                assert (project_dir / "config" / "config.yaml").exists()
                assert (project_dir / "config" / "samples.tsv").exists()
            finally:
                os.chdir(original_cwd)


class TestWorkflowRun:
    """Tests for the workflow run command."""

    def test_workflow_run_help(self):
        """Test that help works for workflow run command."""
        result = runner.invoke(app, ["workflow", "run", "--help"])
        assert result.exit_code == 0
        assert "Run the Snakemake pipeline" in result.stdout
        assert "--cores" in result.stdout
        assert "--use-docker" in result.stdout

    def test_workflow_run_missing_snakefile(self):
        """Test that error is raised when Snakefile is missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            config_dir = Path(tmpdir) / "config"
            config_dir.mkdir()
            (config_dir / "config.yaml").touch()
            
            result = runner.invoke(
                app,
                ["workflow", "run", "--workflow-dir", tmpdir],
            )
            
            assert result.exit_code == 1
            assert "Snakefile not found" in result.stdout

    def test_workflow_run_missing_config(self):
        """Test that error is raised when config is missing."""
        with tempfile.TemporaryDirectory() as tmpdir:
            result = runner.invoke(
                app,
                ["workflow", "run", "--workflow-dir", tmpdir],
            )
            
            # Should fail with either missing snakefile or config error
            assert result.exit_code == 1
            assert "not found" in result.stdout

    @patch("subprocess.run")
    def test_workflow_run_success(self, mock_run):
        """Test successful workflow run."""
        mock_run.return_value = MagicMock(returncode=0)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create required files
            config_dir = Path(tmpdir) / "config"
            config_dir.mkdir()
            (config_dir / "config.yaml").touch()
            
            workflow_dir = Path(tmpdir) / "workflow"
            workflow_dir.mkdir()
            (workflow_dir / "Snakefile").touch()
            
            result = runner.invoke(
                app,
                ["workflow", "run", "--workflow-dir", tmpdir, "--cores", "4"],
            )
            
            assert result.exit_code == 0
            assert "Running:" in result.stdout
            mock_run.assert_called_once()

    @patch("subprocess.run")
    def test_workflow_run_with_options(self, mock_run):
        """Test workflow run with various options."""
        mock_run.return_value = MagicMock(returncode=0)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create required files
            config_dir = Path(tmpdir) / "config"
            config_dir.mkdir()
            (config_dir / "config.yaml").touch()
            
            workflow_dir = Path(tmpdir) / "workflow"
            workflow_dir.mkdir()
            (workflow_dir / "Snakefile").touch()
            
            result = runner.invoke(
                app,
                [
                    "workflow", "run",
                    "--workflow-dir", tmpdir,
                    "--cores", "8",
                    "--use-docker",
                    "--dry-run",
                    "--verbose",
                ],
            )
            
            assert result.exit_code == 0
            
            # Check that snakemake was called with correct arguments
            call_args = mock_run.call_args[0][0]
            assert "snakemake" in call_args
            assert "--use-docker" in call_args
            assert "-n" in call_args  # dry-run
            assert "-v" in call_args  # verbose
            assert "--cores" in call_args
            assert "8" in call_args

    @patch("subprocess.run", side_effect=FileNotFoundError)
    def test_workflow_run_snakemake_not_found(self, mock_run):
        """Test error when snakemake is not installed."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create required files
            config_dir = Path(tmpdir) / "config"
            config_dir.mkdir()
            (config_dir / "config.yaml").touch()
            
            workflow_dir = Path(tmpdir) / "workflow"
            workflow_dir.mkdir()
            (workflow_dir / "Snakefile").touch()
            
            result = runner.invoke(
                app,
                ["workflow", "run", "--workflow-dir", tmpdir],
            )
            
            assert result.exit_code == 1
            assert "snakemake not found" in result.stdout

    @patch("subprocess.run")
    def test_workflow_run_custom_snakemake_args(self, mock_run):
        """Test passing custom snakemake arguments."""
        mock_run.return_value = MagicMock(returncode=0)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create required files
            config_dir = Path(tmpdir) / "config"
            config_dir.mkdir()
            (config_dir / "config.yaml").touch()
            
            workflow_dir = Path(tmpdir) / "workflow"
            workflow_dir.mkdir()
            (workflow_dir / "Snakefile").touch()
            
            result = runner.invoke(
                app,
                [
                    "workflow", "run",
                    "--workflow-dir", tmpdir,
                    "--snakemake-args", "--forceall --keep-going",
                ],
            )
            
            assert result.exit_code == 0
            
            call_args = mock_run.call_args[0][0]
            assert "--forceall" in call_args
            assert "--keep-going" in call_args


class TestCLIIntegration:
    """Integration tests for the CLI."""

    def test_app_help(self):
        """Test main app help."""
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "scnado" in result.stdout
        assert "find-barcodes" in result.stdout or "find_barcodes" in result.stdout
        assert "fragments" in result.stdout
        assert "workflow" in result.stdout

    def test_workflow_subcommand_help(self):
        """Test workflow subcommand help."""
        result = runner.invoke(app, ["workflow", "--help"])
        assert result.exit_code == 0
        assert "Workflow management" in result.stdout
        assert "init" in result.stdout
        assert "run" in result.stdout

    def test_invalid_command(self):
        """Test that invalid commands fail gracefully."""
        result = runner.invoke(app, ["invalid-command"])
        assert result.exit_code != 0
