"""Optional AI boundary: no remote call or invented interpretation in v0.1."""

import os


def ai_status():
    return {"enabled": False, "reason": "No configured and reviewed interpretation provider; complete analysis remains available.",
            "api_key_present": bool(os.environ.get("EMBRYOOMICS_API_KEY"))}


def aggregate_payload(markers, config):
    """Preview the only data an optional future provider may receive."""
    return {"dataset_accession": config["dataset_accession"], "source_url": config["source_url"],
            "markers": markers.head(50).to_dict(orient="records")}
