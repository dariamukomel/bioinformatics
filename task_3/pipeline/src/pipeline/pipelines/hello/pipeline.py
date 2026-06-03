from kedro.pipeline import Pipeline, node, pipeline

from .nodes import say_hello, save_message


def create_pipeline(**kwargs) -> Pipeline:
    return pipeline(
        [
            node(
                func=say_hello,
                inputs=None,
                outputs="hello_message",
                name="say_hello_node",
            ),
            node(
                func=save_message,
                inputs="hello_message",
                outputs="saved_hello_message",
                name="save_message_node",
            ),
        ]
    )