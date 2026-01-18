from tensorflow.keras.layers import *
import tensorflow as tf

from training_func.custom_model import CustomModel

# --------------------------------------------------

def self_defined_model(dropout_rate, model_name=None):

    '''
        Write your deep neural network model in the USER_model.py file.
        Users must strictly obey the following rules:

        Rules:
            -- (1) The input tensor should match the dimensions provided in the input folders.
            -- (2) The final output tensor shape must be (1, 2).
            -- (3) Proper location of the dropout rate should be placed.
    '''

    # Edit the input layer(s) below

    inputs = Input(shape=(28, 28, 1))

    # Implement your model below

    x = Flatten()(inputs)
    x = Dense(100)(x)
    x = Dropout(dropout_rate)(x)
    
    outputs = Dense(2, activation='softmax')(x)

    # End of implementation

    # If there are miltiple inputs, the 'inputs' parameter below should be editted as expected.
    model = CustomModel(inputs=inputs, outputs=outputs, name=model_name)

    return model
