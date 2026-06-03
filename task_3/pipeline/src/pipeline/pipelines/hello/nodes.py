def say_hello() -> str:
    message = "Hello from Kedro pipeline!"
    print(message)
    return message


def save_message(message: str) -> str:
    print(f"Received message: {message}")
    return message