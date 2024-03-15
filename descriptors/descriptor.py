class Descriptor:
    def calculate(self):
        raise NotImplementedError

    def get_data(self) -> dict:
        raise NotImplementedError