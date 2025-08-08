import requests


class APIClient:
    def __init__(self, base_url: str, polling_interval: int = 1):
        self.base_url = base_url
        self.polling_interval = polling_interval

    def _make_request(self, method: str, endpoint: str, **kwargs) -> requests.Response:
        """Make HTTP request and handle basic error checking."""
        url = f"{self.base_url}/{endpoint}"
        response = requests.request(method, url, **kwargs)
        response.raise_for_status()
        return response
