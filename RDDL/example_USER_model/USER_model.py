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
    # Input shape: (sequence_length, num_channels) for 1D convolution
    # Adjust the shape according to your data dimensions
    
    inputs = Input(shape=(1024, 8))  

    # Implement your model below
    # Architecture: Input -> 1D Conv (ReLU) -> Global Average Pooling -> Dense -> Softmax

    # 1D Convolutional Layer with ReLU activation
    x = Conv1D(filters=64, kernel_size=3, activation='relu', padding='same')(inputs)
    
    # Global Average Pooling Layer
    x = GlobalAveragePooling1D()(x)
    
    # Dropout for regularization
    x = Dropout(dropout_rate)(x)
    
    # Dense Layer
    x = Dense(64, activation='relu')(x)
    
    # Softmax output layer (2 classes)
    outputs = Dense(2, activation='softmax')(x)

    # End of implementation

    # If there are miltiple inputs, the 'inputs' parameter below should be editted as expected.
    model = CustomModel(inputs=inputs, outputs=outputs, name=model_name)

    return model
