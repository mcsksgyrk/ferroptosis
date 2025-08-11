from .base import APIClient


class HPAClient(APIClient):
    def __init__(self):
        super().__init__("https://")
