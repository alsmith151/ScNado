from pathlib import Path
from typing import Optional

import yaml
from pydantic import BaseModel, Field, model_validator
from seqnado.config.configs import GenomeConfig


class CATConfig(BaseModel):
    bin_size: int = 5000
    n_features: int = 50_000


class RNAConfig(BaseModel):
    n_top_genes: int = 2000


class ScnadoConfig(BaseModel):
    samples: Path
    barcode_csv: Optional[Path] = None
    barcode_mismatches: int = 2
    genome: GenomeConfig
    cat: CATConfig = Field(default_factory=CATConfig)
    rna: RNAConfig = Field(default_factory=RNAConfig)
    enable_cat: bool = True
    enable_rna: bool = True
    enable_integration: Optional[bool] = None
    scnado_container: str = "oras://ghcr.io/alsmith151/scnado:latest"

    @model_validator(mode="after")
    def set_integration_default(self) -> "ScnadoConfig":
        if self.enable_integration is None:
            self.enable_integration = self.enable_cat and self.enable_rna
        return self

    @classmethod
    def from_yaml(cls, path: Path) -> "ScnadoConfig":
        with open(path) as f:
            return cls(**yaml.safe_load(f))
