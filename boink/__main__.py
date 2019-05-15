from boink.cli import BoinkRunner
from boink.cdbg import cDBGRunner


def main():
    runner = BoinkRunner()
    runner.add_command('cdbg', cDBGRunner)

    return runner.run()
