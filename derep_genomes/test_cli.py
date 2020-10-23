import pytest
from derep_genomes.__main__ import main


def test_main_template():
    with pytest.raises(SystemExit) as pytest_wrapped_e:
        main()
    assert pytest_wrapped_e.type == SystemExit
